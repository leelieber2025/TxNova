//! Residual splice census: cohort-recurrent CIGAR `N` junctions that are
//! not annotated-gene introns. Harvest ignores treat/control labels.
//! `cohort=silent` is control-absent; `cohort=shared` is also in control;
//! `cohort=cohort` is used when the sheet has no control samples.

use crate::bam;
use crate::coverage::{cigar_introns, inferred_strand, pass_read};
use crate::error::{CoreError, Result};
use crate::gtf::{parse_gtf, write_tsv, GeneBody, Transcript};
use crate::interval::{gap, overlaps};
use crate::sys_mem;
use rust_htslib::bam::record::Record;
use rust_htslib::bam::Read;
use serde::Deserialize;
use std::collections::{HashMap, HashSet};
use std::path::Path;

const MAX_PENDING: usize = 500_000;
const MIN_SUPPORT: u64 = 2;
const MAX_CONTROL_SUPPORT: u64 = 0;

#[derive(Debug, Deserialize)]
struct SampleIn {
    sample_id: String,
    bam: String,
    #[serde(default)]
    group: String,
    #[serde(default)]
    dup_flag_seen: bool,
}

#[derive(Debug, Deserialize)]
struct SamplesJson {
    samples: Vec<SampleIn>,
}

#[derive(Debug, Deserialize)]
struct LeakCfg {
    #[serde(default = "d_mapq")]
    min_mapq: u8,
    #[serde(default)]
    require_unique_nh: bool,
    #[serde(default = "d_rf")]
    strandedness: String,
    #[serde(default = "d_skip")]
    skip_duplicate: String,
    #[serde(default = "d_auto")]
    library_layout: String,
    #[serde(default)]
    threads: usize,
    #[serde(default = "d_reps")]
    min_samples: usize,
    #[serde(default = "d_min_support")]
    min_support: u64,
    #[serde(default = "d_max_ctrl")]
    max_control_support: u64,
}

fn d_mapq() -> u8 {
    10
}
fn d_rf() -> String {
    "rf".into()
}
fn d_skip() -> String {
    "auto".into()
}
fn d_auto() -> String {
    "auto".into()
}
fn d_reps() -> usize {
    2
}
fn d_min_support() -> u64 {
    MIN_SUPPORT
}
fn d_max_ctrl() -> u64 {
    MAX_CONTROL_SUPPORT
}

#[derive(Clone, Debug, Hash, Eq, PartialEq, Ord, PartialOrd)]
struct JuncKey {
    chrom: String,
    start: u64,
    end: u64,
    strand: char,
}

#[derive(Clone, Debug)]
struct MergedHit {
    locus_id: String,
    class: String,
}

fn nh(rec: &Record) -> Option<i32> {
    rec.aux(b"NH").ok().and_then(|a| match a {
        rust_htslib::bam::record::Aux::U8(v) => Some(v as i32),
        rust_htslib::bam::record::Aux::I8(v) => Some(v as i32),
        rust_htslib::bam::record::Aux::U16(v) => Some(v as i32),
        rust_htslib::bam::record::Aux::I16(v) => Some(v as i32),
        rust_htslib::bam::record::Aux::U32(v) => Some(v as i32),
        rust_htslib::bam::record::Aux::I32(v) => Some(v),
        _ => None,
    })
}

fn skip_dup(mode: &str, seen: bool) -> bool {
    match mode {
        "always" => true,
        "never" => false,
        _ => seen,
    }
}

fn mate_still_ahead(rec: &Record, tid: i32, pos: i64) -> bool {
    rec.mtid() > tid || (rec.mtid() == tid && rec.mpos() >= pos)
}

fn tid_names(header: &rust_htslib::bam::HeaderView) -> HashMap<i32, String> {
    (0..header.target_count())
        .map(|tid| {
            (
                tid as i32,
                String::from_utf8_lossy(header.tid2name(tid)).into_owned(),
            )
        })
        .collect()
}

fn add_introns(
    rec: &Record,
    names: &HashMap<i32, String>,
    strandedness: &str,
    counts: &mut HashMap<JuncKey, u64>,
) {
    let Some(chrom) = names.get(&rec.tid()) else {
        return;
    };
    let chrom = chrom.clone();
    let strand = inferred_strand(rec, strandedness).unwrap_or('.');
    let mut seen = HashSet::new();
    for (s, e) in cigar_introns(rec) {
        if !seen.insert((s, e)) {
            continue;
        }
        *counts
            .entry(JuncKey {
                chrom: chrom.clone(),
                start: s,
                end: e,
                strand,
            })
            .or_insert(0) += 1;
    }
}

fn census_one_sample(sample: &SampleIn, cfg: &LeakCfg) -> Result<HashMap<JuncKey, u64>> {
    let mut reader = bam::open_sequential(Path::new(&sample.bam))?;
    let names = tid_names(reader.header());
    let skip = skip_dup(&cfg.skip_duplicate, sample.dup_flag_seen);
    let mut pending: HashMap<Vec<u8>, Record> = HashMap::new();
    let mut counts: HashMap<JuncKey, u64> = HashMap::new();
    let mut last_tid = i32::MIN;

    for rec in reader.records() {
        let rec = rec.map_err(|e| CoreError::fail(format!("BAM read: {e}")))?;
        if !pass_read(&rec, cfg.min_mapq, skip) {
            continue;
        }
        if cfg.require_unique_nh && nh(&rec).unwrap_or(0) != 1 {
            continue;
        }
        let paired = rec.is_paired();
        if paired && cfg.library_layout != "single" {
            if !rec.is_proper_pair() {
                continue;
            }
            let qn = rec.qname().to_vec();
            if let Some(mate) = pending.remove(&qn) {
                let mut frag: HashMap<JuncKey, u64> = HashMap::new();
                add_introns(&rec, &names, &cfg.strandedness, &mut frag);
                add_introns(&mate, &names, &cfg.strandedness, &mut frag);
                for k in frag.into_keys() {
                    *counts.entry(k).or_insert(0) += 1;
                }
            } else if rec.tid() == rec.mtid() && rec.mpos() < rec.pos() {
            } else {
                if rec.tid() != last_tid {
                    pending.retain(|_, stored| stored.mtid() >= rec.tid());
                    last_tid = rec.tid();
                }
                if pending.len() >= MAX_PENDING {
                    pending.retain(|_, stored| mate_still_ahead(stored, rec.tid(), rec.pos()));
                    if pending.len() >= MAX_PENDING {
                        let drop_n = (pending.len() / 2).max(1);
                        let mut ranked: Vec<(i32, i64, Vec<u8>)> = pending
                            .iter()
                            .map(|(k, rec)| (rec.mtid(), rec.mpos(), k.clone()))
                            .collect();
                        ranked.sort_by(|a, b| (a.0, a.1, &a.2).cmp(&(b.0, b.1, &b.2)));
                        for (_, _, k) in ranked.into_iter().take(drop_n) {
                            pending.remove(&k);
                        }
                    }
                }
                pending.insert(qn, rec);
            }
        } else {
            add_introns(&rec, &names, &cfg.strandedness, &mut counts);
        }
    }
    Ok(counts)
}

fn merged_introns(transcripts: &[Transcript]) -> HashMap<(String, u64, u64, char), Vec<MergedHit>> {
    let mut out: HashMap<(String, u64, u64, char), Vec<MergedHit>> = HashMap::new();
    for t in transcripts {
        // merged.gtf is the annotation. Every intron here is known.
        for (s, e) in crate::coverage::transcript_introns(t) {
            out.entry((t.chrom.clone(), s, e, t.strand))
                .or_default()
                .push(MergedHit {
                    locus_id: t.gene_id.clone(),
                    class: "overlap".into(),
                });
        }
    }
    out
}

fn status_of(hits: &[MergedHit]) -> (&'static str, String, String) {
    if hits.is_empty() {
        return ("unassembled", String::new(), String::new());
    }
    let all_u = hits.iter().all(|h| h.class == "u");
    let status = if all_u {
        "assembled_u"
    } else {
        "assembled_other"
    };
    let locus = hits[0].locus_id.clone();
    let class = if all_u {
        "u".into()
    } else {
        hits.iter()
            .find(|h| h.class != "u")
            .map(|h| h.class.clone())
            .unwrap_or_else(|| hits[0].class.clone())
    };
    (status, locus, class)
}

fn nearest_gene(
    chrom: &str,
    start: u64,
    end: u64,
    strand: char,
    genes: &[GeneBody],
    unstranded: bool,
) -> (bool, String, String, String, String) {
    // overlaps_gene = either strand (class u). nearest_* stays same-strand
    // (or any strand if the library is unstranded) for the 5 kb distance gate.
    let mut any_overlap = false;
    let mut best: Option<(&GeneBody, u64)> = None;
    for g in genes {
        if g.chrom != chrom {
            continue;
        }
        let same = unstranded || strand == '.' || g.strand == strand;
        if overlaps(start, end, g.start, g.end) {
            any_overlap = true;
            if same {
                return (
                    true,
                    g.gene_id.clone(),
                    g.gene_name.clone(),
                    "0".into(),
                    g.strand.to_string(),
                );
            }
            continue;
        }
        if !same {
            continue;
        }
        let d = gap(start, end, g.start, g.end);
        match best {
            None => best = Some((g, d)),
            Some((_, bd)) if d < bd => best = Some((g, d)),
            _ => {}
        }
    }
    match best {
        Some((g, d)) => (
            any_overlap,
            g.gene_id.clone(),
            g.gene_name.clone(),
            d.to_string(),
            g.strand.to_string(),
        ),
        None => (
            any_overlap,
            String::new(),
            String::new(),
            String::new(),
            String::new(),
        ),
    }
}

pub fn leak_scan(
    merged_gtf: &str,
    samples_json: &str,
    out_tsv: &str,
    cfg_json: &str,
) -> Result<usize> {
    let cfg: LeakCfg = serde_json::from_str(cfg_json)?;
    let samples: SamplesJson = serde_json::from_str(samples_json)?;
    if samples.samples.is_empty() {
        write_tsv(Path::new(out_tsv), &["chrom"], &[])?;
        return Ok(0);
    }
    let merged = parse_gtf(Path::new(merged_gtf))?;
    let intron_ix = merged_introns(&merged.transcripts);
    let unstranded = cfg.strandedness == "unstranded";

    let per_sample = sys_mem::run_bam_jobs(cfg.threads, samples.samples.len(), |si| {
        census_one_sample(&samples.samples[si], &cfg)
    })?;

    let mut keys: HashSet<JuncKey> = HashSet::new();
    for m in &per_sample {
        keys.extend(m.keys().cloned());
    }

    let ctrl_idx: Vec<usize> = samples
        .samples
        .iter()
        .enumerate()
        .filter(|(_, s)| s.group == "control")
        .map(|(i, _)| i)
        .collect();
    let treat_idx: Vec<usize> = samples
        .samples
        .iter()
        .enumerate()
        .filter(|(_, s)| s.group == "treat")
        .map(|(i, _)| i)
        .collect();

    let mut rows = Vec::new();
    let mut keys_sorted: Vec<JuncKey> = keys.into_iter().collect();
    keys_sorted.sort();
    for key in keys_sorted {
        let counts: Vec<u64> = per_sample
            .iter()
            .map(|m| *m.get(&key).unwrap_or(&0))
            .collect();
        let support_sum: u64 = counts.iter().sum();
        let n_detected = counts.iter().filter(|&&c| c > 0).count();
        if n_detected < cfg.min_samples || support_sum < cfg.min_support {
            continue;
        }
        let control_max = ctrl_idx.iter().map(|&i| counts[i]).max().unwrap_or(0);
        let control_n = ctrl_idx.iter().filter(|&&i| counts[i] > 0).count();
        let treat_sum: u64 = treat_idx.iter().map(|&i| counts[i]).sum();
        let treat_n = treat_idx.iter().filter(|&&i| counts[i] > 0).count();
        let cohort = if ctrl_idx.is_empty() {
            "cohort"
        } else if control_max <= cfg.max_control_support {
            "silent"
        } else {
            "shared"
        };
        let hits = intron_ix
            .get(&(key.chrom.clone(), key.start, key.end, key.strand))
            .cloned()
            .or_else(|| {
                if unstranded || key.strand == '.' {
                    let mut acc = Vec::new();
                    for st in ['+', '-'] {
                        if let Some(h) = intron_ix.get(&(key.chrom.clone(), key.start, key.end, st))
                        {
                            acc.extend(h.iter().cloned());
                        }
                    }
                    if acc.is_empty() {
                        None
                    } else {
                        Some(acc)
                    }
                } else {
                    None
                }
            })
            .unwrap_or_default();
        let (status, locus, class) = status_of(&hits);
        if status == "assembled_other" {
            continue;
        }
        let (ov, gid, gname, dist, gstrand) = nearest_gene(
            &key.chrom,
            key.start,
            key.end,
            key.strand,
            &merged.genes,
            unstranded,
        );
        let mut row = vec![
            key.chrom.clone(),
            key.start.to_string(),
            key.end.to_string(),
            key.strand.to_string(),
            (key.end.saturating_sub(key.start) + 1).to_string(),
            status.to_string(),
            locus,
            class,
            if ov { "true" } else { "false" }.into(),
            gid,
            gname,
            dist,
            gstrand,
            control_max.to_string(),
            support_sum.to_string(),
            n_detected.to_string(),
            treat_sum.to_string(),
            treat_n.to_string(),
            control_n.to_string(),
            cohort.to_string(),
        ];
        for c in &counts {
            row.push(c.to_string());
        }
        rows.push(row);
    }

    let mut header = vec![
        "chrom".into(),
        "start".into(),
        "end".into(),
        "strand".into(),
        "intron_nt".into(),
        "status".into(),
        "merged_locus".into(),
        "merged_class".into(),
        "overlaps_gene".into(),
        "nearest_gene_id".into(),
        "nearest_gene_name".into(),
        "nearest_distance_bp".into(),
        "nearest_strand".into(),
        "control_max".into(),
        "support_sum".into(),
        "n_detected".into(),
        "treat_sum".into(),
        "treat_n_detected".into(),
        "control_n_detected".into(),
        "cohort".into(),
    ];
    for s in &samples.samples {
        header.push(format!("{}_count", s.sample_id));
    }
    let header_refs: Vec<&str> = header.iter().map(|s| s.as_str()).collect();
    write_tsv(Path::new(out_tsv), &header_refs, &rows)?;
    Ok(rows.len())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn status_unassembled_and_u() {
        assert_eq!(status_of(&[]).0, "unassembled");
        let u = [MergedHit {
            locus_id: "MSTRG.1".into(),
            class: "u".into(),
        }];
        let (st, loc, cl) = status_of(&u);
        assert_eq!(st, "assembled_u");
        assert_eq!(loc, "MSTRG.1");
        assert_eq!(cl, "u");
        let mixed = [
            MergedHit {
                locus_id: "MSTRG.1".into(),
                class: "u".into(),
            },
            MergedHit {
                locus_id: "ENSG".into(),
                class: "overlap".into(),
            },
        ];
        assert_eq!(status_of(&mixed).0, "assembled_other");
    }

    fn gb(chrom: &str, strand: char, start: u64, end: u64, name: &str) -> GeneBody {
        GeneBody {
            chrom: chrom.into(),
            strand,
            gene_id: name.into(),
            gene_name: name.into(),
            start,
            end,
        }
    }

    #[test]
    fn antisense_gene_body_sets_overlaps_gene() {
        let genes = [
            gb("chr1", '+', 100, 500, "Col"),
            gb("chr1", '-', 10_000, 11_000, "Far"),
        ];
        let (ov, _id, name, dist, gstrand) = nearest_gene("chr1", 200, 300, '-', &genes, false);
        assert!(ov, "opposite-strand gene body is an overlap");
        assert_eq!(name, "Far");
        assert_eq!(gstrand, "-");
        assert_eq!(dist, (10_000 - 300 - 1).to_string());
    }

    #[test]
    fn same_strand_overlap_still_zero_distance() {
        let genes = [gb("chr1", '-', 100, 500, "Hit")];
        let (ov, _id, name, dist, _) = nearest_gene("chr1", 200, 300, '-', &genes, false);
        assert!(ov);
        assert_eq!(name, "Hit");
        assert_eq!(dist, "0");
    }
}
