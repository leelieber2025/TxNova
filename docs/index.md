# TxNova Documentation

TxNova finds **experimental-group-specific novel intergenic transcripts** from
bulk RNA-seq. Input is aligned BAMs. Output is a candidate table and a report.

Treat-recurrent residual splices become locus models. Class, counts, junctions,
and splice-bridging evidence are recomputed from the BAM and the annotation.

## The pipeline, in one picture

```text
 coordinate-sorted, indexed BAMs (STAR / HISAT2)
             │
             ▼
    ┌──────────────────────┐
    │ 1. Known universe      │  slim gene/transcript/exon from the annotation
    └──────────────────────┘
             │
             ▼
    ┌──────────────────────┐
    │ 2. Leak + residual     │  treat-recurrent, control-silent CIGAR N
    │                        │  clustered into RSDL loci; clip off gene bodies
    └──────────────────────┘
             │
             ▼
    ┌──────────────────────┐
    │ 3. Classify            │  universe = annotation + residual
    │    (intergenic, u)     │  u = no exon overlap, no gene body either strand
    └──────────────────────┘
             │
             ▼
    ┌──────────────────────┐
    │ 4. Quantify            │  Rust re-counts every universe locus from BAM
    └──────────────────────┘
             │
             ▼
    ┌──────────────────────┐
    │ 5. Structure + gates   │  splice, distance, valley, TPM, then DE
    └──────────────────────┘
             │  optional
             ▼
    ┌──────────────────────┐
    │ 6. Coding              │  ORF, hexamer LLR, Fickett; fold / orphan
    └──────────────────────┘
```

Per-sample GTFs are stubs. Leak and quantify scan the BAMs.

## Where to go

| Goal | Page |
|------|------|
| Install (`pip install txnova`) | [Installation](installation.md) |
| Run your first analysis end to end | [Quickstart](quickstart.md) |
| Turn raw FASTQs into a TxNova-ready BAM | [Data preparation](data_preparation.md) |
| Every `config.yaml` field, what it does, and its default | [Configuration reference](configuration.md) |
| What's in `output_dir/` and what each column means | [Output reference](outputs.md) |
| Errors, empty tables, and common gotchas | [FAQ / Troubleshooting](faq.md) |
| A worked public-data example (ENCODE BMDM ± Lipid A) | [Public mouse smoke test](PUBLIC_MOUSE.md) |

## The three commands

```bash
txnova init -c config.yaml --samples samples.tsv   # write starter files
txnova preflight -c config.yaml                     # validate BAMs/FASTA/GTF/sample sheet
txnova run -c config.yaml                            # leak → residual → classify → quantify → gate → (DE, coding)
```

`txnova preflight` is cheap and catches most mistakes (mismatched contigs,
un-indexed BAMs, mixed aligners, wrong `group`/`strandedness` values) before
you spend time scanning BAMs. `txnova run` runs preflight again automatically.
There is also `txnova report` to rebuild Markdown/HTML from an existing run.

## What you need before you start

- Coordinate-sorted, **indexed** BAMs from **STAR** or **HISAT2** (not a mix
  of both), at least 1 control and 1 treat sample (≥2 each if you want
  differential expression).
- A genome FASTA (with a `.fai` index) and a matching GTF annotation whose
  contig names agree with the BAM's `@SQ` lines.

TxNova starts from BAMs you already aligned.

## Status

TxNova is `0.1.0`, pre-release, mouse-first (`species: mouse`; the packaged
hexamer table is mouse). Install with `pip install txnova` (see
[Installation](installation.md)). The CLI surface (`init`,
`preflight`, `run`, `report`) and the `config.yaml` schema are the stable,
documented interface; the Python modules under `txnova.*` are internal.
