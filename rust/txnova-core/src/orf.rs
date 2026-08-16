//! Longest ATG…stop ORF on the spliced representative.

use crate::coding::{coding_label, fickett_score, fmt_score, HexamerTable};
use crate::error::{CoreError, Result};
use crate::fasta::FastaIndex;
use crate::gtf::{parse_gtf, write_tsv, Transcript};
use serde::Deserialize;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

const STOPS: &[&[u8]] = &[b"TAG", b"TAA", b"TGA"];

#[derive(Deserialize)]
struct Cfg {
    #[serde(default = "d50")]
    min_orf_aa: usize,
    hexamer_table: String,
    #[serde(default)]
    hexamer_coding_min: f64,
    #[serde(default)]
    hexamer_noncoding_max: f64,
}

fn d50() -> usize {
    50
}

struct OrfHit {
    aa: String,
    nt: Vec<u8>,
    complete: bool,
}

fn splice_seq(t: &Transcript, fa: &FastaIndex) -> Result<Vec<u8>> {
    let mut seq = Vec::new();
    let mut exons = t.exons.clone();
    exons.sort_by_key(|e| e.start);
    for e in &exons {
        seq.extend(fa.fetch(&t.chrom, e.start, e.end)?);
    }
    if t.strand == '-' {
        seq = seq.into_iter().rev().map(rc).collect();
    }
    Ok(seq)
}

fn rc(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' | b'U' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        x => x,
    }
}

fn translate(nt: &[u8]) -> String {
    fn aa(c: &[u8]) -> char {
        match c {
            b"TTT" | b"TTC" => 'F',
            b"TTA" | b"TTG" | b"CTT" | b"CTC" | b"CTA" | b"CTG" => 'L',
            b"ATT" | b"ATC" | b"ATA" => 'I',
            b"ATG" => 'M',
            b"GTT" | b"GTC" | b"GTA" | b"GTG" => 'V',
            b"TCT" | b"TCC" | b"TCA" | b"TCG" | b"AGT" | b"AGC" => 'S',
            b"CCT" | b"CCC" | b"CCA" | b"CCG" => 'P',
            b"ACT" | b"ACC" | b"ACA" | b"ACG" => 'T',
            b"GCT" | b"GCC" | b"GCA" | b"GCG" => 'A',
            b"TAT" | b"TAC" => 'Y',
            b"CAT" | b"CAC" => 'H',
            b"CAA" | b"CAG" => 'Q',
            b"AAT" | b"AAC" => 'N',
            b"AAA" | b"AAG" => 'K',
            b"GAT" | b"GAC" => 'D',
            b"GAA" | b"GAG" => 'E',
            b"TGT" | b"TGC" => 'C',
            b"TGG" => 'W',
            b"CGT" | b"CGC" | b"CGA" | b"CGG" | b"AGA" | b"AGG" => 'R',
            b"GGT" | b"GGC" | b"GGA" | b"GGG" => 'G',
            _ => 'X',
        }
    }
    nt.chunks(3)
        .filter(|c| c.len() == 3 && !STOPS.contains(c))
        .map(aa)
        .collect()
}

fn longest_orf(seq: &[u8], min_aa: usize) -> Option<OrfHit> {
    let mut best: Option<Vec<u8>> = None;
    let mut complete = false;
    for frame in 0..3 {
        let mut i = frame;
        while i + 3 <= seq.len() {
            if &seq[i..i + 3] == b"ATG" {
                let mut j = i + 3;
                let mut found_stop = false;
                while j + 3 <= seq.len() {
                    if STOPS.contains(&&seq[j..j + 3]) {
                        found_stop = true;
                        let orf = seq[i..j].to_vec();
                        if orf.len() / 3 >= min_aa && best.as_ref().map(|b| orf.len() > b.len()).unwrap_or(true)
                        {
                            best = Some(orf);
                            complete = true;
                        }
                        break;
                    }
                    j += 3;
                }
                if !found_stop {
                    let orf = seq[i..seq.len() - (seq.len() - i) % 3].to_vec();
                    if orf.len() / 3 >= min_aa && best.as_ref().map(|b| orf.len() > b.len()).unwrap_or(true)
                    {
                        best = Some(orf);
                        complete = false;
                    }
                }
            }
            i += 3;
        }
    }
    best.map(|nt| {
        let aa = translate(&nt);
        OrfHit { aa, nt, complete }
    })
}

pub fn scan_orfs(
    fasta: &str,
    merged_gtf: &str,
    representatives_tsv: &str,
    out_tsv: &str,
    peptides_fa: &str,
    cfg_json: &str,
) -> Result<usize> {
    let cfg: Cfg = serde_json::from_str(cfg_json)?;
    if cfg.hexamer_table.is_empty() {
        return Err(CoreError::fail("coding.hexamer_table is empty"));
    }
    let _ = cfg.min_orf_aa;
    let hex = HexamerTable::load(Path::new(&cfg.hexamer_table))?;
    let fa = FastaIndex::open(Path::new(fasta))?;
    let parsed = parse_gtf(Path::new(merged_gtf))?;
    let mut by_tid: HashMap<String, &Transcript> = HashMap::new();
    for t in &parsed.transcripts {
        by_tid.insert(t.transcript_id.clone(), t);
    }
    let f = File::open(representatives_tsv)
        .map_err(|e| CoreError::fail(format!("open reps: {e}")))?;
    let mut reps = Vec::new();
    for (i, line) in BufReader::new(f).lines().enumerate() {
        let line = line?;
        if i == 0 && line.contains("transcript_id") {
            continue;
        }
        let mut p = line.split('\t');
        let loc = p.next().unwrap_or("").to_string();
        let tid = p.next().unwrap_or("").to_string();
        if !tid.is_empty() {
            reps.push((loc, tid));
        }
    }
    let mut rows = Vec::new();
    let pep_path = Path::new(peptides_fa);
    if let Some(p) = pep_path.parent() {
        std::fs::create_dir_all(p)?;
    }
    let mut pep = File::create(pep_path)?;
    let mut n_ok = 0usize;
    for (loc, tid) in reps {
        let Some(t) = by_tid.get(&tid) else { continue };
        let seq = splice_seq(t, &fa)?;
        // Report the real longest ATG… ORF. min_orf_aa is the require_orf
        // filter in Python, not a reporting floor — short residual models
        // otherwise show 0 and hide hexamer/Fickett.
        let hit = longest_orf(&seq, 1);
        let (n_aa, complete, score, fickett, label) = match &hit {
            Some(h) => {
                let score = hex.score(&h.nt);
                (
                    h.aa.len(),
                    h.complete,
                    score,
                    Some(fickett_score(&h.nt)),
                    coding_label(score, cfg.hexamer_coding_min, cfg.hexamer_noncoding_max),
                )
            }
            None => (0, false, None, None, "noncoding"),
        };
        rows.push(vec![
            loc.clone(),
            tid.clone(),
            n_aa.to_string(),
            if complete { "true" } else { "false" }.into(),
            fmt_score(score),
            fmt_score(fickett),
            label.into(),
        ]);
        if let Some(h) = hit {
            n_ok += 1;
            writeln!(pep, ">{loc}|{tid}|{n_aa}")?;
            for ch in h.aa.as_bytes().chunks(60) {
                writeln!(pep, "{}", std::str::from_utf8(ch).unwrap())?;
            }
        }
    }
    write_tsv(
        Path::new(out_tsv),
        &[
            "locus_id",
            "transcript_id",
            "longest_orf_aa",
            "orf_complete",
            "coding_score",
            "fickett_score",
            "coding_label",
        ],
        &rows,
    )?;
    Ok(n_ok)
}
