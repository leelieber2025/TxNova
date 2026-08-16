# TxNova Documentation

TxNova finds **experimental-group-specific novel intergenic transcripts** from
bulk RNA-seq — genes that current annotation (e.g. GENCODE) doesn't know about,
that show up in your treated/experimental samples, and that are basically
absent in control. It takes aligned BAMs in, and gives you a candidate table
and a report out.

The output table is a **candidate list**, not a claim of new genes. Novel
models come from treat-recurrent residual splices. TxNova recomputes class,
counts, junctions, and splice-bridging evidence in-process from the BAM and
the annotation.

## The pipeline, in one picture

```text
 coordinate-sorted, indexed BAMs (STAR / HISAT2)
             │
             ▼
    ┌──────────────────────┐
    │ 1. Assemble + merge   │  annotation + residual splice loci
    └──────────────────────┘
             │
             ▼
    ┌──────────────────────┐
    │ 2. Classify            │  in-process: which transcripts are class "u"
    │    (intergenic, u)     │  (no exon overlap, no gene body on either strand)
    └──────────────────────┘
             │
             ▼
    ┌──────────────────────┐
    │ 3. Quantify             │  Rust re-counts every merged transcript/gene
    │    (full universe)      │  from the BAM — TPM/counts
    └──────────────────────┘
             │
             ▼
    ┌──────────────────────┐
    │ 4. Structure scan       │  splice canonicity, coverage discontinuity,
    │                         │  bridging junctions, nearest-gene distance
    └──────────────────────┘
             │
             ▼
    ┌──────────────────────┐
    │ 5. Gates (structure,   │  control-silent + treat-recurrent + gene-like
    │    detection, DE)      │  → candidates.tsv
    └──────────────────────┘
             │  optional
             ▼
    ┌──────────────────────┐
    │ 6. Coding                │  ORF scan, hexamer/Fickett score,
    │    (ORF, fold, orphan)   │  3D structure, conservation, Pfam
    └──────────────────────┘
```

Every step after assembly runs in-process (Rust core + Python orchestration).
There is no external assembler binary.

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
txnova run -c config.yaml                            # assemble → classify → quantify → gate → (DE, coding)
```

`txnova preflight` is cheap and catches most mistakes (mismatched contigs,
un-indexed BAMs, mixed aligners, wrong `group`/`strandedness` values) before
you spend time on assembly. Always run it first — `txnova run` runs it again
automatically, but by then you may have already started building samples.tsv
around a bad assumption.

## What you need before you start

- Coordinate-sorted, **indexed** BAMs from **STAR** or **HISAT2** (not a mix
  of both), at least 1 control and 1 treat sample (≥2 each if you want
  differential expression).
- A genome FASTA (with a `.fai` index) and a matching GTF annotation whose
  contig names agree with the BAM's `@SQ` lines.

No FASTQs, no realignment — TxNova starts from BAMs you already produced.

## Status

TxNova is `0.1.0`, pre-release, mouse-first (`species: mouse`; the packaged
hexamer table is mouse). Install with `pip install txnova` (see
[Installation](installation.md)). The CLI surface (`init`,
`preflight`, `run`, `report`) and the `config.yaml` schema are the stable,
documented interface; the Python modules under `txnova.*` are internal.
