# Quickstart

This walks through one full run: write starter files, point them at your
data, validate, run, and read the report. See [Installation](installation.md)
first if you haven't installed TxNova yet.

## 1. What you need in hand

- Coordinate-sorted, **indexed** BAMs (`.bam` + `.bam.bai`) from STAR or
  HISAT2 — at least 1 control and 1 treat sample, all the **same** aligner
  family and the **same** strandedness.
- A genome FASTA with a `.fai` index, and a GTF annotation whose contig names
  match the BAM's `@SQ` lines exactly (same names, same lengths).

If you only have FASTQs, align them first (`STAR` or `HISAT2`) — TxNova does
not align. See [Data preparation](data_preparation.md) for the full recipe,
including how to figure out `strandedness`.

## 2. Write starter files

```bash
txnova init -c config.yaml --samples samples.tsv
```

This writes a `config.yaml` with every default filled in and a `samples.tsv`
with four placeholder rows. It refuses to overwrite either file if it already
exists.

## 3. Fill in the sample sheet

`samples.tsv` is tab-separated with a fixed set of columns:

| Column | Meaning |
|---|---|
| `sample_id` | Unique ID, `[A-Za-z0-9._-]+` only |
| `bam` | Path to the coordinate-sorted, indexed BAM (absolute, or relative to the sheet's directory) |
| `group` | `control` or `treat` — exactly those two strings, no `wt`/`ko`/`case` aliases |
| `strandedness` | `unstranded`, `fr`, or `rf` — must be the same for every sample in the sheet |
| `replicate` | Optional integer. Auto-numbered per group (1, 2, …) if you leave it blank |

```tsv
sample_id	bam	group	strandedness	replicate
ctrl_1	/data/ctrl_1.bam	control	rf	1
ctrl_2	/data/ctrl_2.bam	control	rf	2
treat_1	/data/treat_1.bam	treat	rf	1
treat_2	/data/treat_2.bam	treat	rf	2
```

You need `≥1` control and `≥1` treat sample to run at all, and `≥2` of each
if `de.enabled: true` (the default) so PyDESeq2 has replicates to work with.

## 4. Point the config at your genome

Edit the `genome:` block in `config.yaml`:

```yaml
genome:
  fasta: /path/to/genome.fa            # needs a .fai next to it
  annotation: /path/to/annotation.gtf  # comprehensive, not a "lean" subset
  annotation_source: GENCODE
  annotation_version: M39
  assembly: GRCm39
  naming_annotation: null              # optional; see Configuration reference
```

Use the **comprehensive** annotation your samples were aligned against (or a
matching one), not a trimmed reference-package GTF — a thin annotation makes
real genes look intergenic. See [Configuration
reference](configuration.md#genome) for why `naming_annotation` is separate
from `annotation`.

Everything else in `config.yaml` (quantify/filters/de/coding) has a working
default — leave it alone for a first run.

## 5. Validate before you commit to a run

```bash
txnova preflight -c config.yaml
```

This checks, in order: BAMs exist and are indexed; contig names/lengths agree
across every BAM, the FASTA, and the GTF; all samples share one aligner
family and one library layout; `group` and `strandedness` values are valid;
you have enough samples for the analysis you asked for. It writes
`output_dir/preflight.json` and exits non-zero with a specific message on the
first failure it hits — fix that one thing and re-run.

## 6. Run

```bash
txnova run -c config.yaml
```

```text
run complete: 3 candidates under /path/to/txnova_out
```

A run re-does preflight, writes the annotation-plus-residual universe,
classifies, quantifies every locus, scans structure, applies gates, runs DE
(if enabled), and — if `coding.enabled: true` — scores ORFs and optional
fold/orphan steps. Unchanged stages are skipped (`output_dir/stamps/`);
`--force` ignores stamps and redos everything.

A run against real data can take anywhere from minutes to hours depending on
BAM size and sample count — BAM scans (leak / quantify) and the
network-backed structure/conservation steps dominate. See
[Installation § network
access](installation.md#optional-network-access-for-structure-and-conservation)
if you want to skip the latter.

## 7. Read the output

```text
output_dir/
├── candidates/candidates.tsv    # the final candidate table
├── report/report.md             # same content as a readable Markdown doc
├── report/report.html           # ...and as a self-contained HTML page
└── quantify/                    # full-universe counts / TPM / DE
```

Open `report/report.html` first — it has the run parameters, a funnel table
(how many transcripts survived each stage), and the final candidate table
with everything you need to sanity-check a row: residual locus, junction/bridge
support, TPM in control vs. treat, and DE statistics. Full column reference:
[Output reference](outputs.md).

### If `candidates.tsv` is empty

Real experiments often have zero or a handful of intergenic, treat-specific
loci. To see which stage cut the count:

1. Check `report/report.md`'s **Funnel** section to see which stage cut the
   count to zero (structure gates, control/treat TPM gates, or DE).
2. Look at `candidates/candidates.gates.tsv` — this is the pre-DE view, so if
   it has rows but `candidates.tsv` doesn't, DE is the filter to look at.
3. See the [FAQ](faq.md#candidatestsv-is-empty) for the specific gate
   thresholds and how to loosen them.

## 8. Rebuild just the report

If you've re-inspected outputs by hand, or changed nothing but want fresh
Markdown/HTML without rerunning the pipeline:

```bash
txnova report -c config.yaml
```

## What's next

| Goal | Page |
|------|------|
| Every config field, its default, and what changing it does | [Configuration reference](configuration.md) |
| Full output directory and column reference | [Output reference](outputs.md) |
| A worked example on public ENCODE mouse data | [Public mouse smoke test](PUBLIC_MOUSE.md) |
| Errors and troubleshooting | [FAQ](faq.md) |
