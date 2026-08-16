# Preparing input BAMs

TxNova starts from **already-aligned, coordinate-sorted, indexed BAMs** — it
does not align FASTQs itself. This page is the recipe for getting from raw
reads to a BAM that passes `txnova preflight`, plus how to fill in
`strandedness` correctly (the single most common `samples.tsv` mistake).

Already have BAMs from STAR or HISAT2? Skip to [4. Sort, index, and
verify](#4-sort-index-and-verify).

## 1. Pick an aligner

TxNova accepts **STAR** or **HISAT2** output — but not a mix of the two
within one run (preflight rejects mixed aligner families across samples).

| Aligner | Version used in TxNova's own testing | When |
|---|---|---|
| **STAR** | 2.7.2a | Standard choice; splice-aware, well-suited to a comprehensive GTF (`--sjdbGTFfile`) |
| **HISAT2** | any recent 2.x | Lower memory footprint; also fine |

Pick one and use it for every sample in the sheet.

## 2. Build a genome index (once per reference)

**STAR:**

```bash
STAR --runMode genomeGenerate \
     --genomeDir star_index \
     --genomeFastaFiles genome.fa \
     --sjdbGTFfile annotation.gtf \
     --sjdbOverhang 99 \
     --runThreadN 8
```

Use the **same** `genome.fa` / `annotation.gtf` (or a name-and-length
compatible pair) you'll later pass to TxNova's `genome.fasta` /
`genome.annotation` — preflight requires the BAM's `@SQ` contig names and
lengths to match the FASTA `.fai` exactly.

**HISAT2:**

```bash
hisat2-build genome.fa hisat2_index
```

## 3. Align each sample

**STAR** (paired-end example; drop the second `--readFilesIn` argument for
single-end):

```bash
STAR --runMode alignReads \
     --genomeDir star_index \
     --readFilesIn sample_R1.fastq.gz sample_R2.fastq.gz \
     --readFilesCommand zcat \
     --outSAMtype BAM SortedByCoordinate \
     --runThreadN 8 \
     --outFileNamePrefix sample_
```

`--outSAMtype BAM SortedByCoordinate` does the coordinate sort for you in
one pass — no separate `samtools sort` needed.

**HISAT2** (piped straight into a coordinate-sorted BAM):

```bash
hisat2 -x hisat2_index -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz -p 8 \
  | samtools sort -@ 8 -o sample_sorted.bam -
```

## 4. Sort, index, and verify

If your BAM isn't already coordinate-sorted:

```bash
samtools sort -@ 8 -o sample_sorted.bam sample.bam
```

Every BAM needs a `.bai` index:

```bash
samtools index sample_sorted.bam
```

Quick sanity check that the contig set matches your FASTA before you write
`samples.tsv`:

```bash
samtools view -H sample_sorted.bam | grep '^@SQ' | head -3
grep '^>' genome.fa | head -3
```

The sequence *names* need to match literally (`chr1` vs `1`, or a GENCODE
vs. Ensembl naming scheme, will both fail preflight with a clear error
naming the mismatching contig).

## 5. Figure out `strandedness`

`samples.tsv` needs one of `unstranded`, `fr`, or `rf` — **the same value
for every sample in the sheet**. This depends on your library prep kit, not
on anything TxNova can infer from the BAM.

| Kit / protocol | `strandedness` |
|---|---|
| Illumina TruSeq Stranded (dUTP-based); most modern "stranded" kits | `rf` (read 2 matches transcript strand) |
| Ligation-based stranded kits (older, less common) | `fr` (read 1 matches transcript strand) |
| Unstranded (e.g. TruSeq RNA, no strand-specific step) | `unstranded` |

If you don't know the kit, check empirically with
[RSeQC](http://rseqc.sourceforge.net/)'s `infer_experiment.py` against your
GTF:

```bash
infer_experiment.py -r annotation.bed12 -i sample_sorted.bam
```

- Mostly `"1++,1--,2+-,2-+"` (or `"++,--"` single-end) → `fr`
- Mostly `"1+-,1-+,2++,2--"` (or `"+-,-+"` single-end) → `rf`
- Roughly 50/50 → `unstranded`

Guessing wrong here doesn't crash the run — it silently biases coverage and
counts, which shows up as spurious `coverage_discontinuity` calls or wrong
`n_exons`/junction assignment. If your funnel numbers look implausible (e.g.
almost nothing survives, or the strand-sensitive structure gates behave
oddly), re-check this first.

## 6. Mark duplicates (optional, but recommended)

`quantify.skip_duplicate: auto` (the default) only skips duplicate-flagged
(`0x400`) reads if it sees the flag actually used in the BAM — preflight
warns per sample if it doesn't. If you want duplicate-aware counting, mark
duplicates before handing the BAM to TxNova:

```bash
samtools markdup -@ 8 sample_sorted.bam sample_markdup.bam
samtools index sample_markdup.bam
```

This isn't required — TxNova runs fine either way — but on PCR-heavy
libraries it changes `control_max_tpm` / `treat_median_tpm` gate outcomes.

## 7. Write the sample sheet and validate

```tsv
sample_id	bam	group	strandedness	replicate
ctrl_1	/data/ctrl_1_sorted.bam	control	rf	1
ctrl_2	/data/ctrl_2_sorted.bam	control	rf	2
treat_1	/data/treat_1_sorted.bam	treat	rf	1
treat_2	/data/treat_2_sorted.bam	treat	rf	2
```

```bash
txnova preflight -c config.yaml
```

Preflight is where all of the above gets checked mechanically: index
present, contig names/lengths agree with the FASTA and across samples, one
aligner family, one library layout, valid `group`/`strandedness` values,
and enough samples for the analysis you configured. See
[Quickstart](quickstart.md) for what happens next, and the
[FAQ](faq.md) for specific preflight error messages.
