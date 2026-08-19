use crate::error::{CoreError, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read as IoRead, Seek, SeekFrom};
use std::path::Path;

/// Parse an existing `.fai`. Never create one.
pub fn load_fai(fasta: &Path) -> Result<Vec<(String, u64)>> {
    let fai = Path::new(&format!("{}.fai", fasta.display())).to_path_buf();
    if !fai.exists() {
        return Err(CoreError::fail(format!(
            "missing FASTA index {} — please run `samtools faidx {}` on that FASTA \
             (TxNova will not write a .fai next to a shared genome)",
            fai.display(),
            fasta.display()
        )));
    }
    let file = File::open(&fai)
        .map_err(|e| CoreError::fail(format!("cannot read {}: {e}", fai.display())))?;
    let mut out = Vec::new();
    for (i, line) in BufReader::new(file).lines().enumerate() {
        let line = line?;
        if line.is_empty() {
            continue;
        }
        let mut parts = line.split('\t');
        let name = parts.next().ok_or_else(|| {
            CoreError::fail(format!("{} line {}: missing contig name", fai.display(), i + 1))
        })?;
        let len: u64 = parts
            .next()
            .ok_or_else(|| {
                CoreError::fail(format!("{} line {}: missing length", fai.display(), i + 1))
            })?
            .parse()
            .map_err(|_| {
                CoreError::fail(format!("{} line {}: invalid length", fai.display(), i + 1))
            })?;
        out.push((name.to_string(), len));
    }
    if out.is_empty() {
        return Err(CoreError::fail(format!("{} is empty", fai.display())));
    }
    let _ = fai; // path used
    Ok(out)
}

pub fn fai_map(fai: &[(String, u64)]) -> HashMap<String, u64> {
    fai.iter().cloned().collect()
}

#[derive(Clone)]
struct FaiRec {
    length: u64,
    offset: u64,
    linebases: u64,
    linewidth: u64,
}

pub struct FastaIndex {
    path: std::path::PathBuf,
    recs: HashMap<String, FaiRec>,
}

impl FastaIndex {
    pub fn open(fasta: &Path) -> Result<Self> {
        let fai = Path::new(&format!("{}.fai", fasta.display())).to_path_buf();
        if !fai.exists() {
            return Err(CoreError::fail(format!(
                "missing FASTA index {} — please run `samtools faidx {}`",
                fai.display(),
                fasta.display()
            )));
        }
        let file = File::open(&fai)
            .map_err(|e| CoreError::fail(format!("cannot read {}: {e}", fai.display())))?;
        let mut recs = HashMap::new();
        for line in BufReader::new(file).lines() {
            let line = line?;
            if line.is_empty() {
                continue;
            }
            let p: Vec<&str> = line.split('\t').collect();
            if p.len() < 5 {
                continue;
            }
            recs.insert(
                p[0].to_string(),
                FaiRec {
                    length: p[1].parse().unwrap_or(0),
                    offset: p[2].parse().unwrap_or(0),
                    linebases: p[3].parse().unwrap_or(1).max(1),
                    linewidth: p[4].parse().unwrap_or(1).max(1),
                },
            );
        }
        Ok(Self {
            path: fasta.to_path_buf(),
            recs,
        })
    }

    /// 1-based inclusive [start, end], genomic strand.
    pub fn fetch(&self, chrom: &str, start: u64, end: u64) -> Result<Vec<u8>> {
        let rec = self.recs.get(chrom).ok_or_else(|| {
            CoreError::fail(format!("FASTA has no contig {chrom}"))
        })?;
        if start == 0 || end < start || end > rec.length {
            return Err(CoreError::fail(format!(
                "bad FASTA interval {chrom}:{start}-{end} (len {})",
                rec.length
            )));
        }
        let start_idx = start - 1;
        let end_idx = end - 1;
        let start_off = rec.offset
            + (start_idx / rec.linebases) * rec.linewidth
            + (start_idx % rec.linebases);
        let end_off = rec.offset
            + (end_idx / rec.linebases) * rec.linewidth
            + (end_idx % rec.linebases);
        let nbytes = (end_off - start_off + 1) as usize;
        let mut f = File::open(&self.path)?;
        f.seek(SeekFrom::Start(start_off))?;
        let mut buf = vec![0u8; nbytes];
        f.read_exact(&mut buf)?;
        let mut out = Vec::with_capacity((end - start + 1) as usize);
        for &b in &buf {
            let c = b.to_ascii_uppercase();
            if c == b'\n' || c == b'\r' {
                continue;
            }
            out.push(c);
        }
        if out.len() != (end - start + 1) as usize {
            return Err(CoreError::fail(format!(
                "FASTA fetch {chrom}:{start}-{end} decoded {} bases, expected {}",
                out.len(),
                end - start + 1
            )));
        }
        Ok(out)
    }

    pub fn dinuc_tx(&self, chrom: &str, genomic_start: u64, plus_strand: bool) -> Result<[u8; 2]> {
        let raw = self.fetch(chrom, genomic_start, genomic_start + 1)?;
        if raw.len() < 2 {
            return Err(CoreError::fail("short dinucleotide fetch"));
        }
        if plus_strand {
            Ok([raw[0], raw[1]])
        } else {
            Ok([comp(raw[1]), comp(raw[0])])
        }
    }
}

fn comp(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' | b'U' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        x => x,
    }
}
