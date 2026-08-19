use crate::error::{CoreError, Result};
use crate::gtf::Transcript;
use rust_htslib::bam::record::Record;
use rust_htslib::bam::Read;
use std::path::Path;

pub fn inferred_strand(rec: &Record, strandedness: &str) -> Option<char> {
    let plus = !rec.is_reverse();
    match strandedness {
        "unstranded" => None,
        "rf" => {
            if rec.is_paired() {
                if rec.is_last_in_template() {
                    Some(if plus { '+' } else { '-' })
                } else {
                    Some(if plus { '-' } else { '+' })
                }
            } else {
                Some(if plus { '-' } else { '+' })
            }
        }
        "fr" => {
            if rec.is_paired() {
                if rec.is_first_in_template() {
                    Some(if plus { '+' } else { '-' })
                } else {
                    Some(if plus { '-' } else { '+' })
                }
            } else {
                Some(if plus { '+' } else { '-' })
            }
        }
        _ => None,
    }
}

pub fn strand_ok(rec: &Record, strandedness: &str, locus_strand: char) -> bool {
    match inferred_strand(rec, strandedness) {
        None => true,
        Some(s) => s == locus_strand || locus_strand == '.',
    }
}

/// 1-based inclusive aligned blocks (M/=/X only).
pub fn aligned_blocks(rec: &Record) -> Vec<(u64, u64)> {
    use rust_htslib::bam::record::Cigar;
    let mut pos = rec.pos() + 1; // 1-based
    if pos < 1 {
        return Vec::new();
    }
    let mut out = Vec::new();
    for c in rec.cigar().iter() {
        match c {
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                let n = *n as i64;
                out.push((pos as u64, (pos + n - 1) as u64));
                pos += n;
            }
            Cigar::Del(n) | Cigar::RefSkip(n) => pos += *n as i64,
            Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
    out
}

pub fn pass_read(rec: &Record, min_mapq: u8, skip_dup: bool) -> bool {
    if rec.is_unmapped() || rec.is_secondary() || rec.is_supplementary() || rec.is_quality_check_failed() {
        return false;
    }
    if skip_dup && rec.is_duplicate() {
        return false;
    }
    rec.mapq() >= min_mapq
}

pub struct ValleyFeat {
    pub locus_mean_depth: f64,
    pub valley_mean: Option<f64>,
    pub valley_possible: Option<bool>,
    pub gap_mean: Option<f64>,
    pub n_dup_flag_seen: u64,
    pub junction_support: Vec<u64>,
    pub bridge_read_count: u64,
}

/// 1-based inclusive introns from CIGAR `N` (not `D`).
pub fn cigar_introns(rec: &Record) -> Vec<(u64, u64)> {
    use rust_htslib::bam::record::Cigar;
    let mut pos = rec.pos() + 1;
    if pos < 1 {
        return Vec::new();
    }
    let mut out = Vec::new();
    for c in rec.cigar().iter() {
        match c {
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => pos += *n as i64,
            Cigar::Del(n) => pos += *n as i64,
            Cigar::RefSkip(n) => {
                let n = *n as i64;
                let start = pos as u64;
                let end = (pos + n - 1) as u64;
                if n >= 1 {
                    out.push((start, end));
                }
                pos += n;
            }
            Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
    out
}

/// 1-based inclusive genomic intron spans. Shared with `splice_features`
/// so `donors` / `junction_support` stay aligned (合同 K).
///
/// Length-1 introns (`e == s`) are kept. Overlapping or abutting exons
/// produce no intron for that pair (`n_introns != n_exons - 1`).
pub fn transcript_introns(t: &Transcript) -> Vec<(u64, u64)> {
    let mut exons = t.exons.clone();
    exons.sort_by_key(|e| e.start);
    let mut out = Vec::new();
    for w in exons.windows(2) {
        let s = w[0].end + 1;
        let e = w[1].start.saturating_sub(1);
        if e >= s {
            out.push((s, e));
        }
    }
    out
}

/// True if a CIGAR `N` joins the novel terminal site facing `gap` to an
/// exon boundary of the nearest gene (合同 J).
///
/// Novel left of the gene: donor is the novel 3′ end (`N` starts at `gap.0`);
/// acceptor is any nearest-gene exon start. Novel right: the mirror.
fn is_bridge(
    ns: &[(u64, u64)],
    t: &Transcript,
    gap: (u64, u64),
    gene_exons: &[(u64, u64)],
) -> bool {
    let (g0, g1) = gap;
    if t.end().saturating_add(1) == g0 {
        for &(es, _) in gene_exons {
            let acc = es.saturating_sub(1);
            if ns.iter().any(|&(s, e)| s == g0 && e == acc) {
                return true;
            }
        }
    } else if t.start() == g1.saturating_add(1) {
        for &(_, ee) in gene_exons {
            let don = ee.saturating_add(1);
            if ns.iter().any(|&(s, e)| s == don && e == g1) {
                return true;
            }
        }
    }
    false
}

fn unique_key(rec: &Record) -> Option<Vec<u8>> {
    match nh(rec) {
        Some(1) => Some(rec.qname().to_vec()),
        Some(_) => None,
        None => {
            let mut k = rec.qname().to_vec();
            k.push(0);
            k.extend(rec.pos().to_le_bytes());
            k.push(0);
            k.extend(format!("{}", rec.cigar()).into_bytes());
            Some(k)
        }
    }
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

pub fn scan_locus_sample(
    reader: &mut rust_htslib::bam::IndexedReader,
    bam_path: &Path,
    chrom: &str,
    t: &Transcript,
    gap: Option<(u64, u64)>,
    gene_exons: &[(u64, u64)],
    valley_bp: u64,
    window_bp: u64,
    strandedness: &str,
    min_mapq: u8,
    skip_dup: bool,
) -> Result<ValleyFeat> {
    let tid = reader
        .header()
        .tid(chrom.as_bytes())
        .ok_or_else(|| CoreError::fail(format!("BAM has no contig {chrom}")))?;

    let fetch_start = match gap {
        Some((g0, _)) => g0.min(t.start()),
        None => t.start(),
    };
    let fetch_end = match gap {
        Some((_, g1)) => g1.max(t.end()),
        None => t.end(),
    };
    // rust-htslib fetch is 0-based half-open
    reader
        .fetch((tid, fetch_start.saturating_sub(1), fetch_end))
        .map_err(|e| CoreError::fail(format!("fetch {}: {e}", bam_path.display())))?;

    let origin = fetch_start;
    let len = (fetch_end - fetch_start + 1) as usize;
    let mut depth = vec![0u32; len];
    let mut n_dup = 0u64;
    let introns = transcript_introns(t);
    let mut junc_seen: Vec<std::collections::HashSet<Vec<u8>>> =
        (0..introns.len()).map(|_| std::collections::HashSet::new()).collect();
    let mut bridge_seen: std::collections::HashSet<Vec<u8>> = std::collections::HashSet::new();

    for rec in reader.records() {
        let rec = rec.map_err(|e| CoreError::fail(format!("BAM read: {e}")))?;
        if rec.is_unmapped() || rec.is_secondary() || rec.is_supplementary() || rec.is_quality_check_failed() {
            continue;
        }
        if rec.mapq() < min_mapq {
            continue;
        }
        if rec.is_duplicate() {
            n_dup += 1;
        }
        if !pass_read(&rec, min_mapq, skip_dup) {
            continue;
        }
        if !strand_ok(&rec, strandedness, t.strand) {
            continue;
        }
        if let Some(key) = unique_key(&rec) {
            let ns = cigar_introns(&rec);
            for (i, iv) in introns.iter().enumerate() {
                if ns.iter().any(|j| j == iv) {
                    junc_seen[i].insert(key.clone());
                }
            }
            if let Some(g) = gap {
                if is_bridge(&ns, t, g, gene_exons) {
                    bridge_seen.insert(key);
                }
            }
        }
        for (s, e) in aligned_blocks(&rec) {
            let a = s.max(fetch_start);
            let b = e.min(fetch_end);
            if a > b {
                continue;
            }
            for p in a..=b {
                let i = (p - origin) as usize;
                depth[i] = depth[i].saturating_add(1);
            }
        }
    }

    let mut exon_sum = 0f64;
    let mut exon_n = 0usize;
    for e in &t.exons {
        let a = e.start.max(fetch_start);
        let b = e.end.min(fetch_end);
        if a > b {
            continue;
        }
        for p in a..=b {
            let i = (p - origin) as usize;
            exon_sum += depth[i] as f64;
            exon_n += 1;
        }
    }
    let locus_mean = if exon_n == 0 { 0.0 } else { exon_sum / exon_n as f64 };

    let (valley_possible, valley_mean) = match gap {
        None => (None, None),
        Some((g0, g1)) if g1 < g0 || (g1 - g0 + 1) < valley_bp => (Some(false), None),
        Some((g0, g1)) => {
            let mut best = f64::MAX;
            let mut start = g0;
            while start + valley_bp - 1 <= g1 {
                let end = start + valley_bp - 1;
                let mut s = 0f64;
                let mut n = 0usize;
                for p in start..=end {
                    if p < fetch_start || p > fetch_end {
                        continue;
                    }
                    s += depth[(p - origin) as usize] as f64;
                    n += 1;
                }
                if n > 0 {
                    best = best.min(s / n as f64);
                }
                start += window_bp.max(1);
            }
            if best == f64::MAX {
                (Some(false), None)
            } else {
                (Some(true), Some(best))
            }
        }
    };

    // Whole-gap mean: a 200 bp empty window next to a transcribed 3′
    // extension (smoke MSTRG.9252) is not an intergenic desert.
    let gap_mean = match gap {
        Some((g0, g1)) if g1 >= g0 => {
            let mut s = 0f64;
            let mut n = 0usize;
            for p in g0..=g1 {
                if p < fetch_start || p > fetch_end {
                    continue;
                }
                s += depth[(p - origin) as usize] as f64;
                n += 1;
            }
            if n == 0 {
                None
            } else {
                Some(s / n as f64)
            }
        }
        _ => None,
    };

    Ok(ValleyFeat {
        locus_mean_depth: locus_mean,
        valley_mean,
        valley_possible,
        gap_mean,
        n_dup_flag_seen: n_dup,
        junction_support: junc_seen.iter().map(|s| s.len() as u64).collect(),
        bridge_read_count: bridge_seen.len() as u64,
    })
}

#[cfg(test)]
mod tests {
    use super::transcript_introns;
    use crate::gtf::{Exon, Transcript};

    #[test]
    fn transcript_introns_genomic() {
        let t = Transcript {
            chrom: "chr1".into(),
            strand: '+',
            gene_id: "g".into(),
            transcript_id: "t".into(),
            gene_name: "g".into(),
            exons: vec![
                Exon { start: 10, end: 20 },
                Exon { start: 40, end: 50 },
                Exon { start: 80, end: 90 },
            ],
        };
        assert_eq!(transcript_introns(&t), vec![(21, 39), (51, 79)]);
    }

    #[test]
    fn one_base_intron_is_kept() {
        let t = Transcript {
            chrom: "chr1".into(),
            strand: '+',
            gene_id: "g".into(),
            transcript_id: "t".into(),
            gene_name: "g".into(),
            exons: vec![Exon { start: 10, end: 20 }, Exon { start: 22, end: 30 }],
        };
        assert_eq!(transcript_introns(&t), vec![(21, 21)]);
    }

    #[test]
    fn abutting_exons_have_no_intron() {
        let t = Transcript {
            chrom: "chr1".into(),
            strand: '+',
            gene_id: "g".into(),
            transcript_id: "t".into(),
            gene_name: "g".into(),
            exons: vec![Exon { start: 10, end: 20 }, Exon { start: 21, end: 30 }],
        };
        assert!(transcript_introns(&t).is_empty());
    }
}
