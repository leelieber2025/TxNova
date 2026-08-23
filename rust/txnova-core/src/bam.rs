use crate::error::{CoreError, Result};
use crate::types::AlignerFamily;
use rust_htslib::bam::{self, IndexedReader, Read};
use sha2::{Digest, Sha256};
use std::fs::File;
use std::io::{Read as IoRead, Seek, SeekFrom};
use std::path::Path;

const BGZF_EOF_LEN: usize = 28;

const POSTPROCESSORS: &[&str] = &[
    "samtools",
    "picard",
    "markduplicates",
    "gatk",
    "sambamba",
    "bedtools",
];

fn is_bgzf_eof_trailer(tail: &[u8; BGZF_EOF_LEN]) -> bool {
    tail[0..4] == [0x1f, 0x8b, 0x08, 0x04]
        && tail[12] == 0x42
        && tail[13] == 0x43
        && tail[16] == 0x1b
        && tail[17] == 0x00
}

pub fn validate_bam_integrity(bam_path: &Path) -> Result<()> {
    let meta = std::fs::metadata(bam_path)
        .map_err(|e| CoreError::fail(format!("cannot stat BAM {}: {e}", bam_path.display())))?;
    let size = meta.len();
    if size < BGZF_EOF_LEN as u64 {
        return Err(CoreError::fail(format!(
            "BAM {} is only {size} bytes; truncated or not a BGZF BAM",
            bam_path.display()
        )));
    }
    let mut f = File::open(bam_path)
        .map_err(|e| CoreError::fail(format!("cannot open BAM {}: {e}", bam_path.display())))?;
    let mut head = [0u8; 2];
    f.read_exact(&mut head).map_err(|e| {
        CoreError::fail(format!(
            "cannot read BAM header bytes {}: {e}",
            bam_path.display()
        ))
    })?;
    if head != [0x1f, 0x8b] {
        return Err(CoreError::fail(format!(
            "BAM {} does not begin with a BGZF gzip header (got {:02x}{:02x})",
            bam_path.display(),
            head[0],
            head[1]
        )));
    }
    f.seek(SeekFrom::End(-(BGZF_EOF_LEN as i64)))
        .map_err(|e| CoreError::fail(format!("cannot seek BAM {}: {e}", bam_path.display())))?;
    let mut tail = [0u8; BGZF_EOF_LEN];
    f.read_exact(&mut tail).map_err(|e| {
        CoreError::fail(format!(
            "cannot read BAM EOF trailer {}: {e}",
            bam_path.display()
        ))
    })?;
    if !is_bgzf_eof_trailer(&tail) {
        return Err(CoreError::fail(format!(
            "BAM {} is missing the BGZF EOF trailer (size {size} bytes); \
             the file is truncated or still being written. Rebuild the index after the BAM is complete.",
            bam_path.display()
        )));
    }
    Ok(())
}

/// Open an indexed BAM. htslib searches `{bam}.bai`, `{bam}.csi`, and `file.bai`.
/// Each reader uses 1 decompression thread so rayon workers do not oversubscribe.
pub fn open_indexed(bam_path: &Path) -> Result<IndexedReader> {
    let mut reader = IndexedReader::from_path(bam_path).map_err(|e| {
        CoreError::fail(format!(
            "cannot open indexed BAM {} (need {0}.bai / {0}.csi / sibling .bai): {e}",
            bam_path.display()
        ))
    })?;
    let _ = reader.set_threads(1);
    Ok(reader)
}

pub fn header_text(reader: &IndexedReader) -> String {
    let header = reader.header();
    String::from_utf8_lossy(header.as_bytes()).into_owned()
}

pub fn parse_sq(header_text: &str) -> Vec<(String, u64)> {
    let mut out = Vec::new();
    for line in header_text.lines() {
        if !line.starts_with("@SQ") {
            continue;
        }
        let mut sn = None;
        let mut ln = None;
        for field in line.split('\t').skip(1) {
            if let Some(v) = field.strip_prefix("SN:") {
                sn = Some(v.to_string());
            } else if let Some(v) = field.strip_prefix("LN:") {
                ln = v.parse().ok();
            }
        }
        if let (Some(sn), Some(ln)) = (sn, ln) {
            out.push((sn, ln));
        }
    }
    out
}

pub fn sq_sha256(sq: &[(String, u64)]) -> String {
    let mut pairs = sq.to_vec();
    pairs.sort_by(|a, b| a.0.cmp(&b.0));
    let mut h = Sha256::new();
    for (sn, ln) in pairs {
        h.update(sn.as_bytes());
        h.update(b"\t");
        h.update(ln.to_string().as_bytes());
        h.update(b"\n");
    }
    format!("{:x}", h.finalize())
}

pub fn is_coordinate_sorted(header_text: &str) -> bool {
    header_text
        .lines()
        .any(|line| line.starts_with("@HD") && line.split('\t').any(|f| f == "SO:coordinate"))
}

#[derive(Debug)]
pub struct LayoutScan {
    pub paired: bool,
    pub mixed: bool,
    pub dup_seen: bool,
}

fn is_postprocessor(id: &str, pn: &str) -> bool {
    let blob = format!("{id} {pn}").to_ascii_lowercase();
    POSTPROCESSORS.iter().any(|tok| blob.contains(tok))
}

pub fn classify_pg_record(id: &str, pn: &str) -> Option<AlignerFamily> {
    if is_postprocessor(id, pn) {
        return None;
    }
    let blob = format!("{id} {pn}").to_ascii_lowercase();
    if blob.contains("bowtie") {
        return Some(AlignerFamily::Bowtie2);
    }
    if blob.contains("minimap") {
        return Some(AlignerFamily::Minimap2);
    }
    if blob.contains("star") {
        return Some(AlignerFamily::Star);
    }
    if blob.contains("hisat") {
        return Some(AlignerFamily::Hisat2);
    }
    Some(AlignerFamily::Other)
}

pub fn parse_pg_records(header_text: &str) -> Vec<(String, String)> {
    let mut out = Vec::new();
    for line in header_text.lines() {
        if !line.starts_with("@PG") {
            continue;
        }
        let mut id = String::new();
        let mut pn = String::new();
        for field in line.split('\t').skip(1) {
            if let Some(v) = field.strip_prefix("ID:") {
                id = v.to_string();
            } else if let Some(v) = field.strip_prefix("PN:") {
                pn = v.to_string();
            }
        }
        out.push((id, pn));
    }
    out
}

/// Contract P §9–14. Any Bowtie2/minimap2 fails even if STAR is also present.
pub fn infer_aligner_family(header_text: &str) -> Result<AlignerFamily> {
    let recs = parse_pg_records(header_text);
    if recs.is_empty() {
        return Err(CoreError::fail(
            "BAM has no @PG records; cannot confirm a splice-aware aligner (STAR or HISAT2)"
                .to_string(),
        ));
    }
    let mut accepted: Vec<AlignerFamily> = Vec::new();
    let mut others = 0usize;
    for (id, pn) in &recs {
        match classify_pg_record(id, pn) {
            None => {}
            Some(f) if f.forbidden() => {
                return Err(CoreError::fail(format!(
                    "BAM @PG contains {} (ID={id} PN={pn}); only STAR or HISAT2 are accepted. \
                     Bowtie2/minimap2 fail even when STAR is also in @PG.",
                    f.as_str()
                )));
            }
            Some(f) if f.accepted() => {
                if !accepted.contains(&f) {
                    accepted.push(f);
                }
            }
            Some(AlignerFamily::Other) => others += 1,
            Some(_) => {}
        }
    }
    match accepted.len() {
        1 => Ok(accepted[0]),
        0 => Err(CoreError::fail(format!(
            "BAM @PG has no STAR/HISAT2 after ignoring post-processors ({} other tool(s)); \
             cannot confirm splice-aware alignment",
            others
        ))),
        _ => Err(CoreError::fail(
            "BAM @PG contains both STAR and HISAT2; mixed aligner families are not allowed"
                .to_string(),
        )),
    }
}

/// Sequential reader for full-file scans (sort / layout). One htslib thread.
pub fn open_sequential(bam_path: &Path) -> Result<bam::Reader> {
    let mut reader = bam::Reader::from_path(bam_path)
        .map_err(|e| CoreError::fail(format!("cannot open BAM {}: {e}", bam_path.display())))?;
    let _ = reader.set_threads(1);
    Ok(reader)
}

pub fn check_coordinate_order_seq(reader: &mut bam::Reader, limit: usize) -> Result<()> {
    let mut prev: Option<(i32, i64)> = None;
    let mut n = 0usize;
    for rec in reader.records() {
        let rec = rec.map_err(|e| CoreError::fail(format!("BAM read error: {e}")))?;
        if rec.is_unmapped() || rec.is_secondary() || rec.is_supplementary() {
            continue;
        }
        let key = (rec.tid(), rec.pos());
        if let Some(p) = prev {
            if key < p {
                return Err(CoreError::fail(
                    "BAM is not coordinate-sorted: primary mapped (tid, pos) decreased \
                     in the first 10000 records"
                        .to_string(),
                ));
            }
        }
        prev = Some(key);
        n += 1;
        if n >= limit {
            break;
        }
    }
    if n < 2 {
        return Err(CoreError::fail(
            "BAM has fewer than 2 primary mapped reads; cannot confirm coordinate sort \
             (unmapped POS=0 records are not evidence)"
                .to_string(),
        ));
    }
    Ok(())
}

pub fn scan_layout_seq(reader: &mut bam::Reader, limit: usize) -> Result<LayoutScan> {
    let mut n_pe = 0usize;
    let mut n_se = 0usize;
    let mut n = 0usize;
    let mut dup_seen = false;
    for rec in reader.records() {
        let rec = rec.map_err(|e| CoreError::fail(format!("BAM read error: {e}")))?;
        if rec.is_unmapped() || rec.is_secondary() || rec.is_supplementary() {
            continue;
        }
        if rec.is_duplicate() {
            dup_seen = true;
        }
        if rec.is_paired() {
            n_pe += 1;
        } else {
            n_se += 1;
        }
        n += 1;
        if n >= limit {
            break;
        }
    }
    if n == 0 {
        return Err(CoreError::fail(
            "BAM has no primary mapped reads; cannot determine library layout".to_string(),
        ));
    }
    let pe_frac = n_pe as f64 / n as f64;
    let se_frac = n_se as f64 / n as f64;
    let mixed = pe_frac > 0.01 && se_frac > 0.01;
    Ok(LayoutScan {
        paired: pe_frac >= se_frac,
        mixed,
        dup_seen,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pg_ignores_samtools() {
        assert!(classify_pg_record("samtools", "samtools").is_none());
        assert!(classify_pg_record("STAR", "STAR") == Some(AlignerFamily::Star));
        assert!(classify_pg_record("bowtie2", "bowtie2") == Some(AlignerFamily::Bowtie2));
    }

    #[test]
    fn star_then_samtools_is_star() {
        let h = "@HD\tVN:1.6\tSO:coordinate\n\
                 @PG\tID:STAR\tPN:STAR\n\
                 @PG\tID:samtools\tPN:samtools\tPP:STAR\n";
        assert_eq!(infer_aligner_family(h).unwrap(), AlignerFamily::Star);
    }

    #[test]
    fn star_plus_bowtie2_fails() {
        let h = "@PG\tID:STAR\tPN:STAR\n@PG\tID:bowtie2\tPN:bowtie2\n";
        assert!(infer_aligner_family(h)
            .unwrap_err()
            .to_string()
            .contains("Bowtie2"));
    }

    #[test]
    fn bowtie2_alone_fails() {
        let h = "@PG\tID:bowtie2\tPN:bowtie2\n";
        assert!(infer_aligner_family(h).is_err());
    }

    #[test]
    fn only_qualimap_fails() {
        let h = "@PG\tID:Qualimap\tPN:Qualimap\n";
        assert!(infer_aligner_family(h).is_err());
    }
}
