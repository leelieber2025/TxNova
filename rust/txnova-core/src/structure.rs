use crate::bam;
use crate::classcode::nearest_genes;
use crate::coverage::{scan_locus_sample, ValleyFeat};
use crate::error::{CoreError, Result};
use crate::fasta::FastaIndex;
use crate::gtf::{parse_gtf, write_tsv, Transcript};
use crate::splice::splice_features;
use crate::sys_mem;
use serde::Deserialize;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
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
struct StructCfg {
    #[serde(default = "d_rf")]
    strandedness: String,
    #[serde(default = "d_mapq")]
    min_mapq: u8,
    #[serde(default)]
    require_unique_nh: bool,
    #[serde(default = "d_skip")]
    skip_duplicate: String,
    #[serde(default = "d_win")]
    discontinuity_window_bp: u64,
    #[serde(default = "d_val")]
    discontinuity_valley_bp: u64,
    #[serde(default)]
    threads: usize,
}

struct LocusPrep {
    locus_id: String,
    t: Transcript,
    donors: String,
    accs: String,
    frac: f64,
    n_introns: usize,
    struct_err: String,
    has_nearest: bool,
    ns: [String; 4],
    na: [String; 4],
    gap: Option<(u64, u64)>,
    gene_exons: Vec<(u64, u64)>,
}

fn d_rf() -> String {
    "rf".into()
}
fn d_mapq() -> u8 {
    10
}
fn d_skip() -> String {
    "auto".into()
}
fn d_win() -> u64 {
    50
}
fn d_val() -> u64 {
    200
}

fn load_reps(path: &Path) -> Result<Vec<(String, String)>> {
    let f = File::open(path).map_err(|e| {
        CoreError::fail(format!(
            "cannot open representatives {}: {e}",
            path.display()
        ))
    })?;
    let mut out = Vec::new();
    for (i, line) in BufReader::new(f).lines().enumerate() {
        let line = line?;
        if i == 0 && line.contains("transcript_id") {
            continue;
        }
        if line.is_empty() {
            continue;
        }
        let mut p = line.split('\t');
        let loc = p.next().unwrap_or("").to_string();
        let tid = p.next().unwrap_or("").to_string();
        if !loc.is_empty() && !tid.is_empty() {
            out.push((loc, tid));
        }
    }
    Ok(out)
}

fn skip_dup_for(mode: &str, seen: bool) -> bool {
    match mode {
        "always" => true,
        "never" => false,
        _ => seen,
    }
}

fn open_gap(a0: u64, a1: u64, b0: u64, b1: u64) -> Option<(u64, u64)> {
    if a1 < b0 {
        let s = a1 + 1;
        let e = b0 - 1;
        if e >= s {
            return Some((s, e));
        }
    } else if b1 < a0 {
        let s = b1 + 1;
        let e = a0 - 1;
        if e >= s {
            return Some((s, e));
        }
    }
    None
}

pub fn structure_scan(
    merged_gtf: &str,
    ref_gtf: &str,
    _class_tsv: &str,
    representatives_tsv: &str,
    fasta: &str,
    samples_json: &str,
    out_tsv: &str,
    cfg_json: &str,
) -> Result<usize> {
    let cfg: StructCfg = serde_json::from_str(cfg_json)?;
    let samples: SamplesJson = serde_json::from_str(samples_json)?;
    let merged = parse_gtf(Path::new(merged_gtf))?;
    let reference = parse_gtf(Path::new(ref_gtf))?;
    let mut by_tid: HashMap<String, &Transcript> = HashMap::new();
    for t in &merged.transcripts {
        by_tid.insert(t.transcript_id.clone(), t);
    }
    let reps = load_reps(Path::new(representatives_tsv))?;
    let fa = FastaIndex::open(Path::new(fasta))?;
    let uns = cfg.strandedness == "unstranded";

    let mut genes_by_chrom: HashMap<String, Vec<&crate::gtf::GeneBody>> = HashMap::new();
    for g in &reference.genes {
        genes_by_chrom.entry(g.chrom.clone()).or_default().push(g);
    }
    let mut exons_by_gene: HashMap<(String, String), Vec<(u64, u64)>> = HashMap::new();
    for rt in &reference.transcripts {
        let slot = exons_by_gene
            .entry((rt.chrom.clone(), rt.gene_id.clone()))
            .or_default();
        for e in &rt.exons {
            slot.push((e.start, e.end));
        }
    }
    let empty_genes: Vec<&crate::gtf::GeneBody> = Vec::new();

    let fmt_near = |n: &Option<(String, String, u64, char, u64, u64)>| -> [String; 4] {
        match n {
            Some((id, name, d, s, _, _)) => {
                [id.clone(), name.clone(), d.to_string(), s.to_string()]
            }
            None => [String::new(), String::new(), String::new(), String::new()],
        }
    };

    let mut preps: Vec<LocusPrep> = Vec::new();
    for (locus_id, tid) in &reps {
        let t = match by_tid.get(tid) {
            Some(t) => (*t).clone(),
            None => continue,
        };
        let uns_t = uns || t.strand == '.';
        let chrom_genes = genes_by_chrom.get(&t.chrom).unwrap_or(&empty_genes);
        let (near_same, near_any) = nearest_genes(&t, chrom_genes, uns_t);
        let (donors, accs, frac, n_introns, structure_error) = match splice_features(&t, &fa) {
            Ok(v) => v,
            Err(e) => (
                String::new(),
                String::new(),
                0.0,
                crate::coverage::transcript_introns(&t).len(),
                Some(e.to_string()),
            ),
        };
        let struct_err = structure_error.unwrap_or_default();
        let nearest_for_gap = if uns_t {
            near_any.clone()
        } else {
            near_same.clone()
        };
        // has_nearest follows the same gene used for gap (stranded = same-strand only).
        let has_nearest = nearest_for_gap.is_some();
        let gap = nearest_for_gap
            .as_ref()
            .and_then(|(_, _, _, _, gs, ge)| open_gap(t.start(), t.end(), *gs, *ge));
        let mut gene_exons: Vec<(u64, u64)> = nearest_for_gap
            .as_ref()
            .and_then(|(gid, _, _, _, _, _)| exons_by_gene.get(&(t.chrom.clone(), gid.clone())))
            .cloned()
            .unwrap_or_default();
        if gene_exons.is_empty() {
            if let Some((_, _, _, _, gs, ge)) = nearest_for_gap.as_ref() {
                gene_exons.push((*gs, *ge));
            }
        }
        preps.push(LocusPrep {
            locus_id: locus_id.clone(),
            t,
            donors,
            accs,
            frac,
            n_introns,
            struct_err,
            has_nearest,
            ns: fmt_near(&near_same),
            na: fmt_near(&near_any),
            gap,
            gene_exons,
        });
    }

    let mut rows = Vec::new();
    if samples.samples.is_empty() {
        for p in &preps {
            rows.push(empty_sample_row(
                &p.locus_id,
                &p.t.transcript_id,
                &p.t,
                &p.donors,
                &p.accs,
                p.frac,
                p.n_introns,
                p.has_nearest,
                &p.ns,
                &p.na,
                &p.struct_err,
            ));
        }
    } else {
        let n_s = samples.samples.len();
        let mut slots: Vec<Vec<String>> = vec![Vec::new(); preps.len() * n_s];
        let scanned = sys_mem::run_bam_jobs(cfg.threads, n_s, |si| {
            scan_one_sample(si, &samples.samples[si], &preps, &cfg).map(|(_, rows)| rows)
        })?;
        for (si, sample_rows) in scanned.into_iter().enumerate() {
            for (li, row) in sample_rows.into_iter().enumerate() {
                slots[li * n_s + si] = row;
            }
        }
        rows = slots;
    }

    write_tsv(
        Path::new(out_tsv),
        &[
            "locus_id",
            "sample_id",
            "transcript_id",
            "chrom",
            "start",
            "end",
            "strand",
            "n_exons",
            "length_nt",
            "donors",
            "acceptors",
            "canonical_splice_fraction",
            "n_introns",
            "has_nearest",
            "nearest_gene_id",
            "nearest_gene_name",
            "nearest_distance_bp",
            "nearest_strand",
            "nearest_any_gene_id",
            "nearest_any_gene_name",
            "nearest_any_distance_bp",
            "nearest_any_strand",
            "locus_mean_depth",
            "valley_mean",
            "valley_possible",
            "gap_mean_depth",
            "n_dup_flag_seen",
            "junction_support",
            "bridge_read_count",
            "structure_error",
        ],
        &rows,
    )?;
    Ok(reps.len())
}

fn scan_one_sample(
    si: usize,
    s: &SampleIn,
    preps: &[LocusPrep],
    cfg: &StructCfg,
) -> Result<(usize, Vec<Vec<String>>)> {
    let mut sample_rows = Vec::with_capacity(preps.len());
    if s.bam.is_empty() {
        for p in preps {
            sample_rows.push(feature_row(p, &s.sample_id, None));
        }
        return Ok((si, sample_rows));
    }
    let mut reader = bam::open_indexed(Path::new(&s.bam))?;
    let skip = skip_dup_for(&cfg.skip_duplicate, s.dup_flag_seen);
    for p in preps {
        let feat = scan_locus_sample(
            &mut reader,
            Path::new(&s.bam),
            &p.t.chrom,
            &p.t,
            p.gap,
            &p.gene_exons,
            cfg.discontinuity_valley_bp,
            cfg.discontinuity_window_bp,
            &cfg.strandedness,
            cfg.min_mapq,
            skip,
            cfg.require_unique_nh,
        )?;
        sample_rows.push(feature_row(p, &s.sample_id, Some(feat)));
    }
    Ok((si, sample_rows))
}

fn feature_row(p: &LocusPrep, sample_id: &str, feat: Option<ValleyFeat>) -> Vec<String> {
    let (mean, vmean, vposs, gmean, ndup, junc, bridge) = match feat {
        Some(f) => (
            format!("{:.6}", f.locus_mean_depth),
            f.valley_mean
                .map(|v| format!("{:.6}", v))
                .unwrap_or_default(),
            f.valley_possible
                .map(|b| if b { "true" } else { "false" }.to_string())
                .unwrap_or_default(),
            f.gap_mean.map(|v| format!("{:.6}", v)).unwrap_or_default(),
            f.n_dup_flag_seen.to_string(),
            f.junction_support
                .iter()
                .map(|n| n.to_string())
                .collect::<Vec<_>>()
                .join(","),
            f.bridge_read_count.to_string(),
        ),
        None => (
            String::new(),
            String::new(),
            String::new(),
            String::new(),
            "0".into(),
            String::new(),
            "0".into(),
        ),
    };
    // no nearest → empty valley/gap (not false). nearest but no open gap
    // (overlap / abut) → valley_possible=false (合同 S §1 vs §3).
    let (vmean, vposs, gmean) = if !p.has_nearest {
        (String::new(), String::new(), String::new())
    } else if p.gap.is_none() {
        (String::new(), "false".into(), String::new())
    } else {
        (vmean, vposs, gmean)
    };
    vec![
        p.locus_id.clone(),
        sample_id.to_string(),
        p.t.transcript_id.clone(),
        p.t.chrom.clone(),
        p.t.start().to_string(),
        p.t.end().to_string(),
        p.t.strand.to_string(),
        p.t.n_exons().to_string(),
        p.t.spliced_len().to_string(),
        p.donors.clone(),
        p.accs.clone(),
        if p.n_introns == 0 {
            String::new()
        } else {
            format!("{:.6}", p.frac)
        },
        p.n_introns.to_string(),
        if p.has_nearest { "true" } else { "false" }.into(),
        p.ns[0].clone(),
        p.ns[1].clone(),
        p.ns[2].clone(),
        p.ns[3].clone(),
        p.na[0].clone(),
        p.na[1].clone(),
        p.na[2].clone(),
        p.na[3].clone(),
        mean,
        vmean,
        vposs,
        gmean,
        ndup,
        junc,
        bridge,
        p.struct_err.clone(),
    ]
}

fn empty_sample_row(
    locus_id: &str,
    tid: &str,
    t: &Transcript,
    donors: &str,
    accs: &str,
    frac: f64,
    n_introns: usize,
    has_nearest: bool,
    ns: &[String; 4],
    na: &[String; 4],
    structure_error: &str,
) -> Vec<String> {
    vec![
        locus_id.to_string(),
        String::new(),
        tid.to_string(),
        t.chrom.clone(),
        t.start().to_string(),
        t.end().to_string(),
        t.strand.to_string(),
        t.n_exons().to_string(),
        t.spliced_len().to_string(),
        donors.to_string(),
        accs.to_string(),
        if n_introns == 0 {
            String::new()
        } else {
            format!("{frac:.6}")
        },
        n_introns.to_string(),
        if has_nearest { "true" } else { "false" }.into(),
        ns[0].clone(),
        ns[1].clone(),
        ns[2].clone(),
        ns[3].clone(),
        na[0].clone(),
        na[1].clone(),
        na[2].clone(),
        na[3].clone(),
        String::new(),
        String::new(),
        String::new(),
        String::new(),
        "0".into(),
        String::new(),
        "0".into(),
        structure_error.to_string(),
    ]
}
