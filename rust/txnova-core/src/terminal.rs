//! Treat-coverage terminal exons for residual splice loci.
//!
//! Junctions fix the intron. Terminal exons walk outward from each splice
//! site while treat depth stays above a floor (Cufflinks / StringTie style).
//! Gene-body clipping stays in Python.

use crate::bam;
use crate::coverage::{aligned_blocks, pass_read, strand_ok};
use crate::error::{CoreError, Result};
use crate::fasta;
use crate::gtf::write_tsv;
use rust_htslib::bam::Read;
use serde::Deserialize;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

#[derive(Debug, Deserialize)]
struct SampleIn {
    bam: String,
    #[serde(default)]
    group: String,
    strandedness: String,
    #[serde(default)]
    dup_flag_seen: bool,
}

#[derive(Debug, Deserialize)]
struct SamplesJson {
    samples: Vec<SampleIn>,
}

#[derive(Debug, Deserialize)]
struct Cfg {
    #[serde(default = "d_mapq")]
    min_mapq: u8,
    #[serde(default = "d_skip")]
    skip_duplicate: String,
    #[serde(default = "d_max")]
    max_terminal_nt: u64,
    #[serde(default = "d_depth")]
    min_depth: u32,
    #[serde(default = "d_gap")]
    max_gap_nt: u64,
    #[serde(default)]
    fasta: Option<String>,
}

fn d_mapq() -> u8 {
    10
}
fn d_skip() -> String {
    "auto".into()
}
fn d_max() -> u64 {
    2000
}
fn d_depth() -> u32 {
    1
}
fn d_gap() -> u64 {
    20
}

fn skip_dup(mode: &str, seen: bool) -> bool {
    match mode {
        "always" => true,
        "never" => false,
        _ => seen,
    }
}

/// Walk from `origin` along `dir` (±1). `depth[i]` is position `window_start + i`.
pub fn walk_terminal(
    depth: &[u32],
    window_start: u64,
    origin: u64,
    dir: i64,
    max_nt: u64,
    min_depth: u32,
    max_gap: u64,
) -> u64 {
    let mut last = origin;
    let mut gap = 0u64;
    for step in 1..=max_nt {
        let p = if dir > 0 {
            origin.saturating_add(step)
        } else {
            match origin.checked_sub(step) {
                Some(x) if x >= 1 => x,
                _ => break,
            }
        };
        if p < window_start {
            break;
        }
        let idx = (p - window_start) as usize;
        if idx >= depth.len() {
            break;
        }
        if depth[idx] >= min_depth {
            last = p;
            gap = 0;
        } else {
            gap += 1;
            if gap > max_gap {
                break;
            }
        }
    }
    last
}

fn add_depth(
    reader: &mut rust_htslib::bam::IndexedReader,
    bam_path: &Path,
    chrom: &str,
    window_start: u64,
    window_end: u64,
    depth: &mut [u32],
    strandedness: &str,
    locus_strand: char,
    min_mapq: u8,
    skip: bool,
) -> Result<()> {
    let tid = match reader.header().tid(chrom.as_bytes()) {
        Some(t) => t,
        None => return Ok(()),
    };
    reader
        .fetch((tid, window_start.saturating_sub(1), window_end))
        .map_err(|e| CoreError::fail(format!("fetch {}: {e}", bam_path.display())))?;
    for rec in reader.records() {
        let rec = rec.map_err(|e| CoreError::fail(format!("BAM read {}: {e}", bam_path.display())))?;
        if !pass_read(&rec, min_mapq, skip) {
            continue;
        }
        if !strand_ok(&rec, strandedness, locus_strand) {
            continue;
        }
        for (s, e) in aligned_blocks(&rec) {
            let lo = s.max(window_start);
            let hi = e.min(window_end);
            if hi < lo {
                continue;
            }
            let a = (lo - window_start) as usize;
            let b = (hi - window_start) as usize;
            for d in depth.iter_mut().take(b + 1).skip(a) {
                *d = d.saturating_add(1);
            }
        }
    }
    Ok(())
}

struct Locus {
    id: String,
    chrom: String,
    strand: char,
    intron_start: u64,
    intron_end: u64,
}

fn load_loci(path: &Path) -> Result<Vec<Locus>> {
    let f = File::open(path).map_err(|e| CoreError::fail(format!("open {}: {e}", path.display())))?;
    let mut out = Vec::new();
    for (i, line) in BufReader::new(f).lines().enumerate() {
        let line = line?;
        if i == 0 && line.contains("residual_id") {
            continue;
        }
        if line.trim().is_empty() {
            continue;
        }
        let mut p = line.split('\t');
        let id = p.next().unwrap_or("").to_string();
        let chrom = p.next().unwrap_or("").to_string();
        let strand = p.next().unwrap_or(".").chars().next().unwrap_or('.');
        let intron_start: u64 = p
            .next()
            .unwrap_or("0")
            .parse()
            .map_err(|_| CoreError::fail(format!("bad intron_start on {id}")))?;
        let intron_end: u64 = p
            .next()
            .unwrap_or("0")
            .parse()
            .map_err(|_| CoreError::fail(format!("bad intron_end on {id}")))?;
        if id.is_empty() || chrom.is_empty() || intron_start == 0 || intron_end < intron_start {
            continue;
        }
        out.push(Locus {
            id,
            chrom,
            strand,
            intron_start,
            intron_end,
        });
    }
    Ok(out)
}

pub fn extend_terminals(
    loci_tsv: &str,
    out_tsv: &str,
    samples_json: &str,
    cfg_json: &str,
) -> Result<usize> {
    let cfg: Cfg = serde_json::from_str(cfg_json)?;
    let samples: SamplesJson = serde_json::from_str(samples_json)?;
    if samples.samples.is_empty() {
        return Err(CoreError::fail("extend_terminals: no samples"));
    }
    let loci = load_loci(Path::new(loci_tsv))?;
    let mut lefts = vec![0u64; loci.len()];
    let mut rights = vec![0u64; loci.len()];

    let mut readers: Vec<(rust_htslib::bam::IndexedReader, bool, String, String)> = Vec::new();
    for sample in &samples.samples {
        let skip = skip_dup(&cfg.skip_duplicate, sample.dup_flag_seen);
        let reader = bam::open_indexed(Path::new(&sample.bam))?;
        readers.push((reader, skip, sample.strandedness.clone(), sample.bam.clone()));
    }

    let contig_len = match &cfg.fasta {
        Some(p) => fasta::fai_map(&fasta::load_fai(Path::new(p))?),
        None => std::collections::HashMap::new(),
    };

    for (i, loc) in loci.iter().enumerate() {
        let reach = cfg.max_terminal_nt.saturating_add(1);
        let win_s = loc.intron_start.saturating_sub(reach).max(1);
        let mut win_e = loc.intron_end.saturating_add(reach);
        if let Some(&ln) = contig_len
            .get(&loc.chrom)
            .or_else(|| contig_len.get(loc.chrom.strip_prefix("chr").unwrap_or(&loc.chrom)))
        {
            win_e = win_e.min(ln);
        }
        let mut depth = vec![0u32; (win_e - win_s + 1) as usize];
        for (reader, skip, strandedness, bam_path) in readers.iter_mut() {
            add_depth(
                reader,
                Path::new(bam_path),
                &loc.chrom,
                win_s,
                win_e,
                &mut depth,
                strandedness,
                loc.strand,
                cfg.min_mapq,
                *skip,
            )?;
        }
        let left_origin = loc.intron_start.saturating_sub(1).max(1);
        let right_origin = loc.intron_end.saturating_add(1);
        lefts[i] = walk_terminal(
            &depth,
            win_s,
            left_origin,
            -1,
            cfg.max_terminal_nt,
            cfg.min_depth,
            cfg.max_gap_nt,
        );
        rights[i] = walk_terminal(
            &depth,
            win_s,
            right_origin,
            1,
            cfg.max_terminal_nt,
            cfg.min_depth,
            cfg.max_gap_nt,
        );
    }

    let mut rows = Vec::with_capacity(loci.len());
    for (i, loc) in loci.iter().enumerate() {
        rows.push(vec![
            loc.id.clone(),
            loc.chrom.clone(),
            loc.strand.to_string(),
            loc.intron_start.to_string(),
            loc.intron_end.to_string(),
            lefts[i].to_string(),
            rights[i].to_string(),
        ]);
    }
    write_tsv(
        Path::new(out_tsv),
        &[
            "residual_id",
            "chrom",
            "strand",
            "intron_start",
            "intron_end",
            "left_start",
            "right_end",
        ],
        &rows,
    )?;
    Ok(loci.len())
}

#[cfg(test)]
mod tests {
    use super::walk_terminal;

    #[test]
    fn walk_stops_after_gap() {
        // positions 10..20, origin 15, depth high at 15-13 then a 3-gap then a blip
        let start = 10u64;
        let depth = vec![0, 0, 3, 3, 3, 2, 0, 0, 0, 1, 0]; // 10..20
        let last = walk_terminal(&depth, start, 15, -1, 20, 1, 2);
        assert_eq!(last, 12);
        let last_r = walk_terminal(&depth, start, 15, 1, 20, 1, 2);
        assert_eq!(last_r, 15);
    }

    #[test]
    fn walk_extends_while_covered() {
        let start = 1u64;
        let depth = vec![0, 1, 1, 1, 1, 1, 0, 0]; // pos 1..8
        let last = walk_terminal(&depth, start, 2, 1, 20, 1, 0);
        assert_eq!(last, 6);
    }
}
