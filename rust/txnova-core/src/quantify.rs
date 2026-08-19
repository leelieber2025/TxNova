use crate::bam;
use crate::coverage::{aligned_blocks, inferred_strand, pass_read};
use crate::error::{CoreError, Result};
use crate::gtf::{parse_gtf, write_tsv, Transcript};
use crate::interval::{Interval, IvIndex};
use crate::sys_mem;
use rust_htslib::bam::record::Record;
use rust_htslib::bam::Read;
use serde::Deserialize;
use std::collections::{HashMap, HashSet};
use std::path::Path;

#[derive(Debug, Deserialize)]
struct SampleIn {
    sample_id: String,
    bam: String,
    #[serde(default)]
    dup_flag_seen: bool,
}

#[derive(Debug, Deserialize)]
struct SamplesJson {
    samples: Vec<SampleIn>,
}

#[derive(Debug, Clone, Deserialize)]
struct QuantCfg {
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

#[derive(Clone)]
struct ExonHit {
    transcript_id: String,
    gene_id: String,
    strand: char,
}

struct TxMeta {
    gene_id: String,
    length: u64,
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

fn compatible_hits(
    rec: &Record,
    idx: &HashMap<i32, &IvIndex<ExonHit>>,
    tid: i32,
    strandedness: &str,
) -> HashSet<(String, String)> {
    let Some(ix) = idx.get(&tid) else {
        return HashSet::new();
    };
    let inf = inferred_strand(rec, strandedness);
    let mut out = HashSet::new();
    for (s, e) in aligned_blocks(rec) {
        for hit in ix.overlapping(s, e) {
            let ok = match inf {
                None => true,
                Some(st) => st == hit.data.strand || hit.data.strand == '.',
            };
            if ok {
                out.insert((hit.data.gene_id.clone(), hit.data.transcript_id.clone()));
            }
        }
    }
    out
}

fn count_fragment(
    hits: &HashSet<(String, String)>,
    locus_c: &mut HashMap<String, f64>,
    tx_c: &mut HashMap<String, f64>,
) {
    let genes: HashSet<&str> = hits.iter().map(|(g, _)| g.as_str()).collect();
    if genes.len() != 1 {
        return;
    }
    let gid = genes.iter().next().unwrap().to_string();
    *locus_c.entry(gid.clone()).or_insert(0.0) += 1.0;
    let txs: Vec<&str> = hits
        .iter()
        .filter(|(g, _)| g == &gid)
        .map(|(_, t)| t.as_str())
        .collect();
    let k = txs.len().max(1) as f64;
    for t in txs {
        *tx_c.entry(t.to_string()).or_insert(0.0) += 1.0 / k;
    }
}

fn skip_dup(mode: &str, seen: bool) -> bool {
    match mode {
        "always" => true,
        "never" => false,
        _ => seen,
    }
}

fn build_exon_index_by_chrom(transcripts: &[Transcript]) -> HashMap<String, IvIndex<ExonHit>> {
    let mut by_chrom: HashMap<String, Vec<Interval<ExonHit>>> = HashMap::new();
    for t in transcripts {
        for e in &t.exons {
            by_chrom.entry(t.chrom.clone()).or_default().push(Interval {
                start: e.start,
                end: e.end,
                data: ExonHit {
                    transcript_id: t.transcript_id.clone(),
                    gene_id: t.gene_id.clone(),
                    strand: t.strand,
                },
            });
        }
    }
    by_chrom
        .into_iter()
        .map(|(k, v)| (k, IvIndex::from_intervals(v)))
        .collect()
}

fn tid_index<'a>(
    header: &rust_htslib::bam::HeaderView,
    by_chrom: &'a HashMap<String, IvIndex<ExonHit>>,
) -> HashMap<i32, &'a IvIndex<ExonHit>> {
    let mut out = HashMap::new();
    for (chrom, iv) in by_chrom {
        if let Some(tid) = header.tid(chrom.as_bytes()) {
            out.insert(tid as i32, iv);
        }
    }
    out
}

fn mate_still_ahead(rec: &Record, tid: i32, pos: i64) -> bool {
    rec.mtid() > tid || (rec.mtid() == tid && rec.mpos() >= pos)
}

/// Bound for unpaired PE first-mates. STAR BAMs are coordinate-sorted;
/// mates whose coordinate is already behind the scan cursor are evicted.
const MAX_PENDING: usize = 500_000;

fn quantify_one_sample(
    sample: &SampleIn,
    cfg: &QuantCfg,
    by_chrom: &HashMap<String, IvIndex<ExonHit>>,
) -> Result<(HashMap<String, f64>, HashMap<String, f64>, u64)> {
    let mut reader = bam::open_sequential(Path::new(&sample.bam))?;
    let idx = tid_index(reader.header(), by_chrom);
    let skip = skip_dup(&cfg.skip_duplicate, sample.dup_flag_seen);
    let mut pending: HashMap<Vec<u8>, Record> = HashMap::new();
    let mut locus_c: HashMap<String, f64> = HashMap::new();
    let mut tx_c: HashMap<String, f64> = HashMap::new();
    let mut last_tid = i32::MIN;
    let mut n_pending_drop = 0u64;

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
                let mut hits = compatible_hits(&rec, &idx, rec.tid(), &cfg.strandedness);
                hits.extend(compatible_hits(&mate, &idx, mate.tid(), &cfg.strandedness));
                count_fragment(&hits, &mut locus_c, &mut tx_c);
            } else if rec.tid() == rec.mtid() && (rec.mpos() as i64) < rec.pos() {
                // Mate is already behind the cursor and was never stored
                // (failed filters). Do not keep this end forever.
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
                            n_pending_drop += 1;
                        }
                    }
                }
                pending.insert(qn, rec);
            }
        } else {
            let hits = compatible_hits(&rec, &idx, rec.tid(), &cfg.strandedness);
            count_fragment(&hits, &mut locus_c, &mut tx_c);
        }
    }

    Ok((locus_c, tx_c, n_pending_drop))
}

pub fn quantify_gtf(
    merged_gtf: &str,
    samples_json: &str,
    out_dir: &str,
    cfg_json: &str,
) -> Result<(usize, usize, u64)> {
    let cfg: QuantCfg = serde_json::from_str(cfg_json)?;
    let samples: SamplesJson = serde_json::from_str(samples_json)?;
    let parsed = parse_gtf(Path::new(merged_gtf))?;

    let mut tx_meta: HashMap<String, TxMeta> = HashMap::new();
    let mut gene_txs: HashMap<String, Vec<&Transcript>> = HashMap::new();
    for t in &parsed.transcripts {
        tx_meta.insert(
            t.transcript_id.clone(),
            TxMeta {
                gene_id: t.gene_id.clone(),
                length: t.spliced_len(),
            },
        );
        gene_txs.entry(t.gene_id.clone()).or_default().push(t);
    }
    let mut gene_len: HashMap<String, u64> = HashMap::new();
    for (g, ts) in &gene_txs {
        gene_len.insert(g.clone(), Transcript::union_len(ts));
    }

    let sample_ids: Vec<String> = samples.samples.iter().map(|s| s.sample_id.clone()).collect();
    let mut locus_mat: HashMap<String, Vec<f64>> = HashMap::new();
    let mut tx_mat: HashMap<String, Vec<f64>> = HashMap::new();
    for t in &parsed.transcripts {
        tx_mat.insert(t.transcript_id.clone(), vec![0.0; sample_ids.len()]);
        locus_mat
            .entry(t.gene_id.clone())
            .or_insert_with(|| vec![0.0; sample_ids.len()]);
    }

    let mut n_pending_drop_total = 0u64;
    if !samples.samples.is_empty() {
        let by_chrom = build_exon_index_by_chrom(&parsed.transcripts);
        let results = sys_mem::run_bam_jobs(cfg.threads, samples.samples.len(), |si| {
            quantify_one_sample(&samples.samples[si], &cfg, &by_chrom)
        })?;
        for (si, slot) in results.into_iter().enumerate() {
            let (locus_c, tx_c, drops) = slot;
            n_pending_drop_total += drops;
            for (g, c) in locus_c {
                if let Some(row) = locus_mat.get_mut(&g) {
                    row[si] = c;
                }
            }
            for (t, c) in tx_c {
                if let Some(row) = tx_mat.get_mut(&t) {
                    row[si] = c;
                }
            }
        }
    }

    // TPM on full tables
    let mut locus_rpk: Vec<(String, Vec<f64>, Vec<f64>)> = Vec::new();
    for (g, counts) in &locus_mat {
        let ell = *gene_len.get(g).unwrap_or(&1) as f64;
        let ell = ell.max(1.0);
        let rpk: Vec<f64> = counts.iter().map(|c| c / ell).collect();
        locus_rpk.push((g.clone(), counts.clone(), rpk));
    }
    let n_s = sample_ids.len();
    let mut loc_denom = vec![0.0; n_s];
    for (_, _, rpk) in &locus_rpk {
        for i in 0..n_s {
            loc_denom[i] += rpk[i];
        }
    }

    let mut tx_rpk: Vec<(String, Vec<f64>, Vec<f64>)> = Vec::new();
    for (tid, counts) in &tx_mat {
        let ell = tx_meta.get(tid).map(|m| m.length).unwrap_or(1) as f64;
        let ell = ell.max(1.0);
        let rpk: Vec<f64> = counts.iter().map(|c| c / ell).collect();
        tx_rpk.push((tid.clone(), counts.clone(), rpk));
    }
    let mut tx_denom = vec![0.0; n_s];
    for (_, _, rpk) in &tx_rpk {
        for i in 0..n_s {
            tx_denom[i] += rpk[i];
        }
    }

    let out = Path::new(out_dir);
    std::fs::create_dir_all(out)?;

    let mut loc_count_rows = Vec::new();
    let mut loc_tpm_rows = Vec::new();
    locus_rpk.sort_by(|a, b| a.0.cmp(&b.0));
    for (g, counts, rpk) in &locus_rpk {
        let mut cr = vec![g.clone()];
        let mut tr = vec![g.clone()];
        for i in 0..n_s {
            cr.push(format!("{:.0}", counts[i]));
            let tpm = if loc_denom[i] > 0.0 {
                rpk[i] / loc_denom[i] * 1e6
            } else {
                0.0
            };
            tr.push(format!("{tpm:.6}"));
        }
        loc_count_rows.push(cr);
        loc_tpm_rows.push(tr);
    }
    let mut loc_header = vec!["locus_id".to_string()];
    loc_header.extend(sample_ids.clone());
    let loc_h: Vec<&str> = loc_header.iter().map(|s| s.as_str()).collect();
    write_tsv(&out.join("locus_counts.tsv"), &loc_h, &loc_count_rows)?;
    write_tsv(&out.join("locus_tpm.tsv"), &loc_h, &loc_tpm_rows)?;

    let mut tx_count_rows = Vec::new();
    let mut tx_tpm_rows = Vec::new();
    tx_rpk.sort_by(|a, b| a.0.cmp(&b.0));
    for (tid, counts, rpk) in &tx_rpk {
        let gid = tx_meta.get(tid).map(|m| m.gene_id.clone()).unwrap_or_default();
        let mut cr = vec![tid.clone(), gid.clone()];
        let mut tr = vec![tid.clone(), gid];
        for i in 0..n_s {
            cr.push(format!("{:.6}", counts[i]));
            let tpm = if tx_denom[i] > 0.0 {
                rpk[i] / tx_denom[i] * 1e6
            } else {
                0.0
            };
            tr.push(format!("{tpm:.6}"));
        }
        tx_count_rows.push(cr);
        tx_tpm_rows.push(tr);
    }
    let mut tx_header = vec!["transcript_id".to_string(), "locus_id".to_string()];
    tx_header.extend(sample_ids);
    let tx_h: Vec<&str> = tx_header.iter().map(|s| s.as_str()).collect();
    write_tsv(&out.join("transcript_counts.tsv"), &tx_h, &tx_count_rows)?;
    write_tsv(&out.join("transcript_tpm.tsv"), &tx_h, &tx_tpm_rows)?;

    Ok((parsed.transcripts.len(), locus_mat.len(), n_pending_drop_total))
}
