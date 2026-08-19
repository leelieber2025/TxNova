use crate::error::{CoreError, Result};
use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

#[derive(Clone, Debug)]
pub struct Exon {
    pub start: u64,
    pub end: u64,
}

#[derive(Clone, Debug)]
pub struct Transcript {
    pub chrom: String,
    pub strand: char,
    pub gene_id: String,
    pub transcript_id: String,
    pub gene_name: String,
    pub exons: Vec<Exon>,
}

impl Transcript {
    pub fn start(&self) -> u64 {
        self.exons.iter().map(|e| e.start).min().unwrap_or(0)
    }

    pub fn end(&self) -> u64 {
        self.exons.iter().map(|e| e.end).max().unwrap_or(0)
    }

    pub fn n_exons(&self) -> usize {
        self.exons.len()
    }

    pub fn spliced_len(&self) -> u64 {
        self.exons.iter().map(|e| e.end - e.start + 1).sum()
    }

    pub fn union_len(transcripts: &[&Transcript]) -> u64 {
        if transcripts.is_empty() {
            return 0;
        }
        let mut ivs: Vec<(u64, u64)> = transcripts
            .iter()
            .flat_map(|t| t.exons.iter().map(|e| (e.start, e.end)))
            .collect();
        ivs.sort_unstable();
        let mut tot = 0u64;
        let mut cs = ivs[0].0;
        let mut ce = ivs[0].1;
        for &(s, e) in &ivs[1..] {
            if s <= ce + 1 {
                ce = ce.max(e);
            } else {
                tot += ce - cs + 1;
                cs = s;
                ce = e;
            }
        }
        tot + (ce - cs + 1)
    }
}

#[derive(Clone, Debug)]
pub struct GeneBody {
    pub chrom: String,
    pub strand: char,
    pub gene_id: String,
    pub gene_name: String,
    pub start: u64,
    pub end: u64,
}

/// Collect seqnames from gene / transcript / exon rows. Comments skipped.
pub fn gtf_seqnames(gtf: &Path) -> Result<BTreeSet<String>> {
    let parsed = parse_gtf(gtf)?;
    if parsed.seqnames.is_empty() {
        return Err(CoreError::fail(format!(
            "GTF {} has no gene/transcript/exon rows",
            gtf.display()
        )));
    }
    Ok(parsed.seqnames)
}

#[derive(Debug)]
pub struct ParsedGtf {
    pub transcripts: Vec<Transcript>,
    pub genes: Vec<GeneBody>,
    pub seqnames: BTreeSet<String>,
}

pub fn parse_gtf(gtf: &Path) -> Result<ParsedGtf> {
    let file = File::open(gtf)
        .map_err(|e| CoreError::fail(format!("cannot open GTF {}: {e}", gtf.display())))?;
    let mut tx_map: BTreeMap<String, Transcript> = BTreeMap::new();
    let mut gene_feat: BTreeMap<(String, String), GeneBody> = BTreeMap::new();
    let mut exon_span: BTreeMap<(String, String), (char, String, u64, u64)> = BTreeMap::new();
    let mut seqnames = BTreeSet::new();

    for (i, line) in BufReader::new(file).lines().enumerate() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 9 {
            return Err(CoreError::fail(format!(
                "{} line {}: expected ≥9 tab-separated columns, got {}",
                gtf.display(),
                i + 1,
                cols.len()
            )));
        }
        let chrom = cols[0];
        let feature = cols[2];
        if !matches!(feature, "gene" | "transcript" | "exon") {
            continue;
        }
        seqnames.insert(chrom.to_string());
        let start: u64 = cols[3].parse().map_err(|_| {
            CoreError::fail(format!("{} line {}: bad GTF start", gtf.display(), i + 1))
        })?;
        let end: u64 = cols[4].parse().map_err(|_| {
            CoreError::fail(format!("{} line {}: bad GTF end", gtf.display(), i + 1))
        })?;
        if start > end {
            return Err(CoreError::fail(format!(
                "{} line {}: start {start} > end {end}",
                gtf.display(),
                i + 1
            )));
        }
        let strand = cols[6].chars().next().unwrap_or('.');
        let attrs = cols[8];
        let gene_id = attr(attrs, "gene_id").unwrap_or_default();
        if gene_id.is_empty() {
            return Err(CoreError::fail(format!(
                "{} line {}: gene/transcript/exon row has no gene_id",
                gtf.display(),
                i + 1
            )));
        }
        let gene_name = attr(attrs, "gene_name").unwrap_or_else(|| gene_id.clone());
        let transcript_id = attr(attrs, "transcript_id");

        match feature {
            "gene" => {
                gene_feat.insert(
                    (chrom.to_string(), gene_id.clone()),
                    GeneBody {
                        chrom: chrom.to_string(),
                        strand,
                        gene_id: gene_id.clone(),
                        gene_name,
                        start,
                        end,
                    },
                );
            }
            "transcript" => {
                if let Some(tid) = transcript_id {
                    tx_map.entry(tid.clone()).or_insert(Transcript {
                        chrom: chrom.to_string(),
                        strand,
                        gene_id: gene_id.clone(),
                        transcript_id: tid,
                        gene_name,
                        exons: Vec::new(),
                    });
                }
            }
            "exon" => {
                let tid = match transcript_id {
                    Some(t) if !t.is_empty() => t,
                    _ => format!("{gene_id}.exon"),
                };
                let e = tx_map.entry(tid.clone()).or_insert(Transcript {
                    chrom: chrom.to_string(),
                    strand,
                    gene_id: gene_id.clone(),
                    transcript_id: tid.clone(),
                    gene_name: gene_name.clone(),
                    exons: Vec::new(),
                });
                e.exons.push(Exon { start, end });
                let span = exon_span
                    .entry((chrom.to_string(), gene_id.clone()))
                    .or_insert((strand, gene_name, start, end));
                span.2 = span.2.min(start);
                span.3 = span.3.max(end);
            }
            _ => {}
        }
    }

    let mut transcripts: Vec<Transcript> = tx_map.into_values().filter(|t| !t.exons.is_empty()).collect();
    for t in &mut transcripts {
        t.exons.sort_by_key(|e| e.start);
    }

    let mut genes = Vec::new();
    let mut seen = BTreeSet::new();
    for (key, g) in gene_feat {
        seen.insert(key);
        genes.push(g);
    }
    for ((chrom, gid), (strand, name, start, end)) in exon_span {
        if seen.contains(&(chrom.clone(), gid.clone())) {
            continue;
        }
        genes.push(GeneBody {
            chrom,
            strand,
            gene_id: gid,
            gene_name: name,
            start,
            end,
        });
    }

    Ok(ParsedGtf {
        transcripts,
        genes,
        seqnames,
    })
}

fn attr(attrs: &str, key: &str) -> Option<String> {
    for part in attrs.split(';') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        let mut it = part.splitn(2, char::is_whitespace);
        let k = it.next()?;
        if k != key {
            continue;
        }
        let mut v = it.next().unwrap_or("").trim();
        if v.starts_with('"') && v.ends_with('"') && v.len() >= 2 {
            v = &v[1..v.len() - 1];
        }
        return Some(v.to_string());
    }
    None
}

pub fn write_tsv(path: &Path, header: &[&str], rows: &[Vec<String>]) -> Result<()> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let mut f = File::create(path)
        .map_err(|e| CoreError::fail(format!("cannot write {}: {e}", path.display())))?;
    writeln!(f, "{}", header.join("\t"))?;
    for row in rows {
        writeln!(f, "{}", row.join("\t"))?;
    }
    Ok(())
}

fn strip_chr_prefix(name: &str) -> String {
    let lower = name.to_ascii_lowercase();
    if let Some(rest) = lower.strip_prefix("chr") {
        rest.to_string()
    } else {
        name.to_ascii_lowercase()
    }
}

/// Contract P §8. Returns Ok(n_gtf_seqnames) or a fail message.
pub fn check_gtf_vs_bam(gtf_names: &BTreeSet<String>, bam_names: &BTreeSet<String>) -> Result<usize> {
    if gtf_names.is_empty() {
        return Err(CoreError::fail("GTF has no usable seqnames".to_string()));
    }
    let inter: BTreeSet<_> = gtf_names.intersection(bam_names).cloned().collect();
    if inter.is_empty() {
        let gtf_ex: Vec<_> = gtf_names.iter().take(5).cloned().collect();
        let bam_ex: Vec<_> = bam_names.iter().take(5).cloned().collect();
        return Err(CoreError::fail(format!(
            "GTF seqnames and BAM @SQ have empty intersection \
             (GTF e.g. {gtf_ex:?}, BAM e.g. {bam_ex:?}). \
             Names must match literally; TxNova will not rewrite chr prefixes."
        )));
    }

    let gtf_stripped: BTreeSet<String> = gtf_names.iter().map(|s| strip_chr_prefix(s)).collect();
    let bam_stripped: BTreeSet<String> = bam_names.iter().map(|s| strip_chr_prefix(s)).collect();
    if gtf_stripped == bam_stripped && gtf_names != bam_names {
        let gtf_ex: Vec<_> = gtf_names.iter().take(5).cloned().collect();
        let bam_ex: Vec<_> = bam_names.iter().take(5).cloned().collect();
        return Err(CoreError::fail(format!(
            "GTF vs BAM contig names differ only by a chr/CHR prefix \
             (GTF e.g. {gtf_ex:?}, BAM e.g. {bam_ex:?}). \
             This would classify every transcript as intergenic `u`. Refusing to continue."
        )));
    }

    let extra: Vec<_> = gtf_names.difference(bam_names).take(8).cloned().collect();
    if !extra.is_empty() {
        return Err(CoreError::fail(format!(
            "GTF has seqnames absent from BAM @SQ (e.g. {extra:?}). \
             Use a primary_assembly GTF whose seqnames are a subset of the BAM, \
             e.g. gencode.vM39.primary_assembly.annotation.gtf. Names are not rewritten."
        )));
    }
    Ok(gtf_names.len())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn prefix_conflict() {
        let gtf = BTreeSet::from(["1".into(), "2".into()]);
        let bam = BTreeSet::from(["chr1".into(), "chr2".into()]);
        assert!(check_gtf_vs_bam(&gtf, &bam).is_err());
    }

    #[test]
    fn gtf_extra_contig_fails() {
        let gtf = BTreeSet::from(["chr1".into(), "chr2".into()]);
        let bam = BTreeSet::from(["chr1".into()]);
        assert!(check_gtf_vs_bam(&gtf, &bam)
            .unwrap_err()
            .to_string()
            .contains("absent from BAM"));
    }

    #[test]
    fn subset_gtf_ok() {
        let gtf = BTreeSet::from(["chr1".into()]);
        let bam = BTreeSet::from(["chr1".into(), "chr2".into()]);
        assert_eq!(check_gtf_vs_bam(&gtf, &bam).unwrap(), 1);
    }

    fn write_gtf(name: &str, body: &str) -> std::path::PathBuf {
        let p = std::env::temp_dir().join(name);
        std::fs::write(&p, body).unwrap();
        p
    }

    #[test]
    fn short_line_fails() {
        let p = write_gtf("txnova_short.gtf", "chr1\tX\texon\t1\t10\n");
        assert!(parse_gtf(&p).unwrap_err().to_string().contains("≥9"));
        let _ = std::fs::remove_file(&p);
    }

    #[test]
    fn missing_gene_id_fails() {
        let p = write_gtf(
            "txnova_nogene.gtf",
            "chr1\tX\texon\t1\t10\t.\t+\t.\ttranscript_id \"T1\";\n",
        );
        assert!(parse_gtf(&p).unwrap_err().to_string().contains("gene_id"));
        let _ = std::fs::remove_file(&p);
    }

    #[test]
    fn inverted_coords_fail() {
        let p = write_gtf(
            "txnova_inv.gtf",
            "chr1\tX\texon\t20\t10\t.\t+\t.\tgene_id \"G\"; transcript_id \"T\";\n",
        );
        assert!(parse_gtf(&p).unwrap_err().to_string().contains("start 20 > end 10"));
        let _ = std::fs::remove_file(&p);
    }

    #[test]
    fn same_gene_id_two_chroms_are_separate_bodies() {
        let p = write_gtf(
            "txnova_par.gtf",
            "chrX\tX\texon\t1\t10\t.\t+\t.\tgene_id \"G\"; transcript_id \"TX\";\n\
             chrY\tX\texon\t100\t200\t.\t+\t.\tgene_id \"G\"; transcript_id \"TY\";\n",
        );
        let parsed = parse_gtf(&p).unwrap();
        assert_eq!(parsed.genes.len(), 2);
        let mut bodies: Vec<_> = parsed.genes.iter().map(|g| (g.chrom.as_str(), g.start, g.end)).collect();
        bodies.sort();
        assert_eq!(bodies, vec![("chrX", 1, 10), ("chrY", 100, 200)]);
        let _ = std::fs::remove_file(&p);
    }

    #[test]
    fn exon_only_gtf_builds_gene_body() {
        let p = write_gtf(
            "txnova_exononly.gtf",
            "chr1\tX\texon\t100\t200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
             chr1\tX\texon\t300\t400\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n",
        );
        let parsed = parse_gtf(&p).unwrap();
        assert_eq!(parsed.genes.len(), 1);
        assert_eq!(parsed.genes[0].start, 100);
        assert_eq!(parsed.genes[0].end, 400);
        let _ = std::fs::remove_file(&p);
    }
}
