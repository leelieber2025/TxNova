//! 合同 C — class code.
//!
//! `u(T)` iff no exon–exon overlap on either strand (any strand if unstranded)
//! AND no exon-union vs any gene body on either strand.
//! Prefer missing a dubious overlap over putting overlap into `u`.

use crate::error::Result;
use crate::gtf::{parse_gtf, write_tsv, GeneBody, Transcript};
use crate::interval::{gap, Interval, IvIndex};
use serde::Deserialize;
use std::collections::HashMap;
use std::path::Path;

#[derive(Debug, Deserialize)]
struct ClassifyCfg {
    #[serde(default = "default_strand")]
    strandedness: String,
}

fn default_strand() -> String {
    "rf".into()
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Class {
    U,
    I,
    X,
    Extension,
    Overlap,
}

impl Class {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::U => "u",
            Self::I => "i",
            Self::X => "x",
            Self::Extension => "extension",
            Self::Overlap => "overlap",
        }
    }
}

struct RefIdx {
    exons: HashMap<String, IvIndex<ExonHit>>,
    bodies: HashMap<String, IvIndex<BodyHit>>,
    transcripts: HashMap<String, IvIndex<TxHit>>,
}

#[derive(Clone)]
struct ExonHit {
    strand: char,
    transcript_id: String,
}

#[derive(Clone)]
struct BodyHit {
    strand: char,
    gene_id: String,
    gene_name: String,
}

#[derive(Clone)]
struct TxHit {
    strand: char,
    transcript_id: String,
    start: u64,
    end: u64,
}

fn same_strand_hit(query: char, ref_s: char, unstranded: bool) -> bool {
    if unstranded {
        return true;
    }
    if ref_s == '.' {
        return true;
    }
    query == ref_s
}

fn antisense_hit(query: char, ref_s: char, unstranded: bool) -> bool {
    if unstranded || query == '.' || ref_s == '.' {
        return false;
    }
    query != ref_s
}

fn unstranded_mode(sample_strandedness: &str, query_strand: char) -> bool {
    sample_strandedness == "unstranded" || query_strand == '.'
}

fn build_ref(genes: &[GeneBody], txs: &[Transcript]) -> RefIdx {
    let mut exons: HashMap<String, Vec<Interval<ExonHit>>> = HashMap::new();
    let mut bodies: HashMap<String, Vec<Interval<BodyHit>>> = HashMap::new();
    let mut transcripts: HashMap<String, Vec<Interval<TxHit>>> = HashMap::new();

    for t in txs {
        let span = Interval {
            start: t.start(),
            end: t.end(),
            data: TxHit {
                strand: t.strand,
                transcript_id: t.transcript_id.clone(),
                start: t.start(),
                end: t.end(),
            },
        };
        transcripts.entry(t.chrom.clone()).or_default().push(span);
        for e in &t.exons {
            exons.entry(t.chrom.clone()).or_default().push(Interval {
                start: e.start,
                end: e.end,
                data: ExonHit {
                    strand: t.strand,
                    transcript_id: t.transcript_id.clone(),
                },
            });
        }
    }
    for g in genes {
        bodies.entry(g.chrom.clone()).or_default().push(Interval {
            start: g.start,
            end: g.end,
            data: BodyHit {
                strand: g.strand,
                gene_id: g.gene_id.clone(),
                gene_name: g.gene_name.clone(),
            },
        });
    }

    RefIdx {
        exons: exons
            .into_iter()
            .map(|(k, v)| (k, IvIndex::from_intervals(v)))
            .collect(),
        bodies: bodies
            .into_iter()
            .map(|(k, v)| (k, IvIndex::from_intervals(v)))
            .collect(),
        transcripts: transcripts
            .into_iter()
            .map(|(k, v)| (k, IvIndex::from_intervals(v)))
            .collect(),
    }
}

pub fn classify_transcript(t: &Transcript, idx: &RefIdx, sample_strandedness: &str) -> Class {
    let uns = unstranded_mode(sample_strandedness, t.strand);
    let empty_ex = IvIndex::from_intervals(Vec::new());
    let empty_bd = IvIndex::from_intervals(Vec::new());
    let empty_tx = IvIndex::from_intervals(Vec::new());
    let ex_ix = idx.exons.get(&t.chrom).unwrap_or(&empty_ex);
    let bd_ix = idx.bodies.get(&t.chrom).unwrap_or(&empty_bd);
    let tx_ix = idx.transcripts.get(&t.chrom).unwrap_or(&empty_tx);

    let mut any_exon = false;
    let mut same_exon = false;
    let mut anti_exon = false;
    let mut overlap_len: HashMap<String, u64> = HashMap::new();

    for e in &t.exons {
        for hit in ex_ix.overlapping(e.start, e.end) {
            any_exon = true;
            if same_strand_hit(t.strand, hit.data.strand, uns) {
                same_exon = true;
                let ov = e.end.min(hit.end) - e.start.max(hit.start) + 1;
                *overlap_len.entry(hit.data.transcript_id.clone()).or_insert(0) += ov;
            } else if antisense_hit(t.strand, hit.data.strand, uns) {
                anti_exon = true;
            }
        }
    }

    let t0 = t.start();
    let t1 = t.end();
    let mut same_body = false;
    let mut anti_body = false;
    let mut any_body_any_strand = false;
    // Contract C: body overlap is exon-union vs gene body, not the transcript span.
    for e in &t.exons {
        for hit in bd_ix.overlapping(e.start, e.end) {
            any_body_any_strand = true;
            if same_strand_hit(t.strand, hit.data.strand, uns) {
                same_body = true;
            } else if antisense_hit(t.strand, hit.data.strand, uns) {
                anti_body = true;
            }
        }
    }

    if !any_exon && !any_body_any_strand {
        return Class::U;
    }

    if same_exon {
        // nearest overlapping ref transcript = max shared exon overlap
        let best = overlap_len.into_iter().max_by_key(|(_, n)| *n);
        if let Some((tid, _)) = best {
            for hit in tx_ix.overlapping(t0, t1) {
                if hit.data.transcript_id == tid
                    && (t0 < hit.data.start || t1 > hit.data.end)
                {
                    return Class::Extension;
                }
            }
        }
        return Class::Overlap;
    }
    if anti_exon {
        return Class::X;
    }
    if same_body {
        return Class::I;
    }
    if anti_body {
        return Class::X;
    }
    // body overlap that didn't classify (shouldn't happen) — not u
    Class::Overlap
}

fn exon_structure(t: &Transcript) -> String {
    let mut exons = t.exons.clone();
    exons.sort_by_key(|e| e.start);
    exons
        .iter()
        .map(|e| format!("{}-{}", e.start, e.end))
        .collect::<Vec<_>>()
        .join(",")
}

pub fn classify_gtfs(merged_gtf: &str, ref_gtf: &str, out_tsv: &str, cfg_json: &str) -> Result<(usize, usize)> {
    let cfg: ClassifyCfg = serde_json::from_str(cfg_json)?;
    let merged = parse_gtf(Path::new(merged_gtf))?;
    let reference = parse_gtf(Path::new(ref_gtf))?;
    let idx = build_ref(&reference.genes, &reference.transcripts);

    let mut rows = Vec::new();
    let mut n_u = 0usize;
    for t in &merged.transcripts {
        let class = classify_transcript(t, &idx, &cfg.strandedness);
        if class == Class::U {
            n_u += 1;
        }
        rows.push(vec![
            t.transcript_id.clone(),
            t.gene_id.clone(),
            t.chrom.clone(),
            t.start().to_string(),
            t.end().to_string(),
            t.strand.to_string(),
            t.n_exons().to_string(),
            t.spliced_len().to_string(),
            exon_structure(t),
            class.as_str().to_string(),
            t.gene_name.clone(),
        ]);
    }
    write_tsv(
        Path::new(out_tsv),
        &[
            "transcript_id",
            "gene_id",
            "chrom",
            "start",
            "end",
            "strand",
            "n_exons",
            "length_nt",
            "exon_structure",
            "class",
            "gene_name",
        ],
        &rows,
    )?;
    Ok((merged.transcripts.len(), n_u))
}

/// (gene_id, gene_name, distance_bp, strand, body_start, body_end)
pub type NearestGene = (String, String, u64, char, u64, u64);

pub fn nearest_genes(
    t: &Transcript,
    genes: &[&GeneBody],
    unstranded: bool,
) -> (Option<NearestGene>, Option<NearestGene>) {
    let t0 = t.start();
    let t1 = t.end();
    let mut same: Option<(u64, &GeneBody)> = None;
    let mut any: Option<(u64, &GeneBody)> = None;
    for g in genes {
        if g.chrom != t.chrom {
            continue;
        }
        let d = gap(t0, t1, g.start, g.end);
        match &any {
            None => any = Some((d, g)),
            Some((bd, _)) if d < *bd => any = Some((d, g)),
            _ => {}
        }
        let same_ok = unstranded || t.strand == '.' || g.strand == '.' || t.strand == g.strand;
        if same_ok {
            match &same {
                None => same = Some((d, g)),
                Some((bd, _)) if d < *bd => same = Some((d, g)),
                _ => {}
            }
        }
    }
    let map = |x: Option<(u64, &GeneBody)>| {
        x.map(|(d, g)| {
            (
                g.gene_id.clone(),
                g.gene_name.clone(),
                d,
                g.strand,
                g.start,
                g.end,
            )
        })
    };
    (map(same), map(any))
}
