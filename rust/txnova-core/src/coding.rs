//! Hexamer log-likelihood and Fickett TESTCODE.
//!
//! `coding_score` is mean log2(P_coding / P_noncoding) over in-frame hexamers
//! of the longest ORF. Fickett is a separate column (Fickett 1982 tables as
//! used by CPAT/CPC2). This is not CPAT's logistic model.

use crate::error::{CoreError, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

pub struct HexamerTable {
    llr: HashMap<[u8; 6], f64>,
}

impl HexamerTable {
    pub fn load(path: &Path) -> Result<Self> {
        let f = File::open(path)
            .map_err(|e| CoreError::fail(format!("hexamer table {}: {e}", path.display())))?;
        let mut llr = HashMap::new();
        for (i, line) in BufReader::new(f).lines().enumerate() {
            let line = line?;
            if line.is_empty() {
                continue;
            }
            let mut parts = line.split('\t');
            let hex = parts.next().unwrap_or("");
            if i == 0 && hex.eq_ignore_ascii_case("hexamer") {
                continue;
            }
            let hex_b = hex.as_bytes();
            if hex_b.len() != 6
                || !hex_b
                    .iter()
                    .all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't'))
            {
                return Err(CoreError::fail(format!(
                    "hexamer table {} line {}: expected 6-mer, got {hex:?}",
                    path.display(),
                    i + 1
                )));
            }
            let pc: f64 = parts
                .next()
                .ok_or_else(|| {
                    CoreError::fail(format!(
                        "hexamer table {} line {}: missing coding frequency",
                        path.display(),
                        i + 1
                    ))
                })?
                .parse()
                .map_err(|_| {
                    CoreError::fail(format!(
                        "hexamer table {} line {}: invalid coding frequency",
                        path.display(),
                        i + 1
                    ))
                })?;
            let pn: f64 = parts
                .next()
                .ok_or_else(|| {
                    CoreError::fail(format!(
                        "hexamer table {} line {}: missing noncoding frequency",
                        path.display(),
                        i + 1
                    ))
                })?
                .parse()
                .map_err(|_| {
                    CoreError::fail(format!(
                        "hexamer table {} line {}: invalid noncoding frequency",
                        path.display(),
                        i + 1
                    ))
                })?;
            if !pc.is_finite() || !pn.is_finite() || pc < 0.0 || pn < 0.0 {
                return Err(CoreError::fail(format!(
                    "hexamer table {} line {}: frequencies must be finite and ≥ 0",
                    path.display(),
                    i + 1
                )));
            }
            // Laplace-style floor so a zero on one side stays a strong
            // (non)coding signal instead of dropping the hexamer.
            const PSEUDO: f64 = 1e-8;
            let mut key = [0u8; 6];
            for (j, &b) in hex_b.iter().enumerate() {
                key[j] = b.to_ascii_uppercase();
            }
            llr.insert(key, ((pc + PSEUDO) / (pn + PSEUDO)).log2());
        }
        if llr.is_empty() {
            return Err(CoreError::fail(format!(
                "hexamer table {} has no usable rows",
                path.display()
            )));
        }
        Ok(Self { llr })
    }

    /// Mean in-frame hexamer LLR. `None` if the ORF has no scored hexamer.
    pub fn score(&self, orf_nt: &[u8]) -> Option<f64> {
        let mut sum = 0.0;
        let mut n = 0usize;
        let mut i = 0;
        while i + 6 <= orf_nt.len() {
            let mut key = [0u8; 6];
            let mut ok = true;
            for j in 0..6 {
                let b = orf_nt[i + j].to_ascii_uppercase();
                if matches!(b, b'A' | b'C' | b'G' | b'T') {
                    key[j] = b;
                } else {
                    ok = false;
                    break;
                }
            }
            if ok {
                if let Some(&v) = self.llr.get(&key) {
                    sum += v;
                    n += 1;
                }
            }
            i += 3;
        }
        if n == 0 {
            None
        } else {
            Some(sum / n as f64)
        }
    }
}

/// Fickett 1982 TESTCODE lookup (same tables as CPAT/CPC2).
const POS_CUT: [f64; 10] = [1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1, 0.0];
const POS_PROB: [[f64; 10]; 4] = [
    [0.94, 0.68, 0.84, 0.93, 0.58, 0.68, 0.45, 0.34, 0.20, 0.22], // A
    [0.80, 0.70, 0.70, 0.81, 0.66, 0.48, 0.51, 0.33, 0.30, 0.23], // C
    [0.90, 0.88, 0.74, 0.64, 0.53, 0.48, 0.27, 0.16, 0.08, 0.08], // G
    [0.97, 0.97, 0.91, 0.68, 0.69, 0.44, 0.54, 0.20, 0.09, 0.09], // T
];
const POS_W: [f64; 4] = [0.26, 0.18, 0.31, 0.33];

const CONT_CUT: [f64; 10] = [0.33, 0.31, 0.29, 0.27, 0.25, 0.23, 0.21, 0.19, 0.17, 0.0];
const CONT_PROB: [[f64; 10]; 4] = [
    [0.28, 0.49, 0.44, 0.55, 0.62, 0.49, 0.67, 0.65, 0.81, 0.21], // A
    [0.82, 0.73, 0.43, 0.59, 0.59, 0.32, 0.69, 0.33, 0.30, 0.23], // C
    [0.40, 0.54, 0.47, 0.64, 0.64, 0.73, 0.41, 0.41, 0.33, 0.29], // G
    [0.28, 0.24, 0.39, 0.40, 0.55, 0.75, 0.56, 0.69, 0.51, 0.58], // T
];
const CONT_W: [f64; 4] = [0.11, 0.12, 0.15, 0.14];

fn lookup(value: f64, cuts: &[f64; 10], probs: &[f64; 10]) -> f64 {
    for i in 0..10 {
        if value >= cuts[i] {
            return probs[i];
        }
    }
    probs[9]
}

fn base_idx(b: u8) -> Option<usize> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' | b'U' => Some(3),
        _ => None,
    }
}

pub fn fickett_score(seq: &[u8]) -> f64 {
    if seq.len() < 2 {
        return 0.0;
    }
    let mut pos = [[0u32; 3]; 4];
    for (i, &b) in seq.iter().enumerate() {
        if let Some(bi) = base_idx(b) {
            pos[bi][i % 3] += 1;
        }
    }
    let n = seq.len() as f64;
    let mut score = 0.0;
    for bi in 0..4 {
        let tot = pos[bi][0] + pos[bi][1] + pos[bi][2];
        let content = tot as f64 / n;
        let mx = pos[bi][0].max(pos[bi][1]).max(pos[bi][2]) as f64;
        let mn = pos[bi][0].min(pos[bi][1]).min(pos[bi][2]) as f64;
        let position = mx / (mn + 1.0);
        score += lookup(position, &POS_CUT, &POS_PROB[bi]) * POS_W[bi];
        score += lookup(content, &CONT_CUT, &CONT_PROB[bi]) * CONT_W[bi];
    }
    score
}

/// Mutual exclusive labels, in this order (strict inequalities).
pub fn coding_label(score: Option<f64>, coding_min: f64, noncoding_max: f64) -> &'static str {
    match score {
        None => "no_orf",
        Some(s) if !s.is_finite() => "no_orf",
        Some(s) if s > coding_min => "coding",
        Some(s) if s < noncoding_max => "noncoding",
        Some(_) => "ambiguous",
    }
}

pub fn fmt_score(v: Option<f64>) -> String {
    match v {
        Some(x) if x.is_finite() => format!("{x:.6}"),
        _ => "NA".into(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::Write;

    #[test]
    fn label_contract() {
        assert_eq!(coding_label(None, 0.0, 0.0), "no_orf");
        assert_eq!(coding_label(Some(0.1), 0.0, 0.0), "coding");
        assert_eq!(coding_label(Some(-0.1), 0.0, 0.0), "noncoding");
        assert_eq!(coding_label(Some(0.0), 0.0, 0.0), "ambiguous");
    }

    #[test]
    fn hexamer_mean_llr() {
        let dir = std::env::temp_dir();
        let p = dir.join("txnova_hexamer_test.tsv");
        let mut f = File::create(&p).unwrap();
        writeln!(f, "hexamer\tcoding\tnoncoding").unwrap();
        writeln!(f, "AAAAAA\t0.8\t0.2").unwrap();
        writeln!(f, "TTTTTT\t0.2\t0.8").unwrap();
        drop(f);
        let tab = HexamerTable::load(&p).unwrap();
        let pos = tab.score(b"AAAAAAAA").unwrap(); // one hexamer AAAAAA
        let expect = ((0.8f64 + 1e-8) / (0.2 + 1e-8)).log2();
        assert!((pos - expect).abs() < 1e-12);
        let neg = tab.score(b"TTTTTTTT").unwrap();
        assert!(neg < 0.0);
        assert!(tab.score(b"ATG").is_none());
        let _ = std::fs::remove_file(&p);
    }

    #[test]
    fn zero_frequency_rows_are_kept() {
        let dir = std::env::temp_dir();
        let p = dir.join("txnova_hexamer_zero.tsv");
        let mut f = File::create(&p).unwrap();
        writeln!(f, "hexamer\tcoding\tnoncoding").unwrap();
        writeln!(f, "AAAAAA\t0.8\t0.0").unwrap();
        writeln!(f, "TTTTTT\t0.0\t0.8").unwrap();
        drop(f);
        let tab = HexamerTable::load(&p).unwrap();
        let pos = tab.score(b"AAAAAAAA").unwrap();
        let neg = tab.score(b"TTTTTTTT").unwrap();
        assert!(
            pos > 10.0,
            "coding-only hexamer should be a large +LLR, got {pos}"
        );
        assert!(
            neg < -10.0,
            "noncoding-only hexamer should be a large -LLR, got {neg}"
        );
        let _ = std::fs::remove_file(&p);
    }

    #[test]
    fn fickett_is_finite() {
        let s = fickett_score(b"ATGGCCATTGATAAGCCCAAGTAG");
        assert!(s.is_finite());
        assert!(s > 0.0);
        assert!(s < 2.0);
    }
}
