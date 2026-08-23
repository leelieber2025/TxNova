# Quickstart

This walks through one full run: write starter files, point them at your
data, validate, run, and read the report. See [Installation](installation.md)
first if you haven't installed TxNova yet.

## 1. What you need in hand

- Coordinate-sorted, **indexed** BAMs (`.bam` + `.bam.bai`) from STAR or
  HISAT2 — at least 2 samples, all the **same** aligner family and the
  **same** strandedness. Control and treat labels are optional.
- A genome FASTA with a `.fai` index, and a GTF annotation whose contig names
  match the BAM's `@SQ` lines exactly (same names, same lengths).

If you only have FASTQs, align them first (`STAR` or `HISAT2`) — TxNova does
not align. See [Preparing input BAMs](tutorials/t_prepare_bams.ipynb) for the
recipe, including how to figure out `strandedness`.

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
| `group` | Optional. `control` or `treat` — exactly those two strings, no `wt`/`ko`/`case` aliases. Omit the column (or leave every cell empty) for a catalog-only run. |
| `strandedness` | `unstranded`, `fr`, or `rf` — must be the same for every sample in the sheet |
| `replicate` | Optional integer. Auto-numbered per group (or across the cohort if there is no group) if you leave it blank |

```text
sample_id	bam	group	strandedness	replicate
ctrl_1	/data/ctrl_1.bam	control	rf	1
ctrl_2	/data/ctrl_2.bam	control	rf	2
treat_1	/data/treat_1.bam	treat	rf	1
treat_2	/data/treat_2.bam	treat	rf	2
```

You need `≥2` samples to harvest residual splices. A control-versus-treat
contrast filter (and DE) runs only when the sheet has both groups; DE
additionally needs `≥2` of each and is skipped otherwise.

## 4. Point the config at your genome

Edit the `genome:` block in `config.yaml`:

```yaml
genome:
  fasta: /path/to/genome.fa            # needs a .fai next to it
  annotation: /path/to/annotation.gtf  # comprehensive GENCODE, not a lean subset
  annotation_source: GENCODE
  annotation_version: M39              # or the human GENCODE version
  assembly: GRCm39                     # GRCh38 for human
  naming_annotation: null              # optional; see Configuration reference
```

Mouse: GENCODE comprehensive M39 / GRCm39. Human: GENCODE comprehensive on
GRCh38. `species: auto` (default) infers from the GTF; or set `species: human`.
Only mouse and human are supported.

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

This checks, in order: `group` and `strandedness` values are valid, and you
have enough samples for the analysis you asked for; BAMs exist and are
indexed, with contig names/lengths agreeing across every BAM, the FASTA, and
the GTF; all samples share one aligner family and one library layout. It writes
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
`--force` ignores stamps and reruns every stage.

A run against real data can take anywhere from minutes to hours depending on
BAM size and sample count — BAM scans (leak / quantify) and the
network-backed structure/conservation steps dominate. See
[Installation § network
access](installation.md#optional-network-access-for-structure-and-conservation)
if you want to skip the latter.

## 7. Read the output

```text
output_dir/
├── candidates/residual.tsv             # residual catalog (always)
├── candidates/candidates.tsv           # structure-pass; contrast screen if both groups
├── candidates/candidates.unnamed.tsv   # structure-pass, also in control (TPM)
├── candidates/candidates.shared.tsv    # structure-pass, splice in both groups
├── report/report.html                  # start here
└── quantify/                           # full-universe counts / TPM / DE
```

Open `report/report.html` first, then `candidates/residual.tsv`. Without
control and treat, `candidates.tsv` **is** the structure-pass catalog.
With both groups, it is the treat-detected / control-silent screen; unnamed
and shared hold structure-pass loci that are also in control. Column
reference: [Output reference](outputs.md#the-three-tables).
Main-task walkthrough: [Residual catalog](tutorials/t_residual_catalog.ipynb).

### If `candidates.tsv` is empty

Without a contrast, an empty `candidates.tsv` means every residual failed
a structure gate (splice, distance, valley, bridge, RepeatMasker). Look at
`candidates/residual.tsv` first.

With a contrast, an empty **screen** is common. Check unnamed (control TPM
still on) and shared (splice in both groups) before loosening gates.

If those are empty too:

1. Check `report/report.md`'s **Funnel** for which stage cut the count
   (structure, control/treat TPM, or DE).
2. `candidates/candidates.gates.tsv` is the pre-DE view of the screen.
   Rows there but not in `candidates.tsv` means DE filtered them.
3. See the [FAQ](faq.md#candidatestsv-is-empty).

## 8. Rebuild just the report

If you've re-inspected outputs by hand, or changed nothing but want fresh
Markdown/HTML without rerunning the pipeline:

```bash
txnova report -c config.yaml
```

## What's next

| Goal | Page |
|------|------|
| FASTQ → BAM, strandedness, preflight | [Preparing input BAMs](tutorials/t_prepare_bams.ipynb) |
| Every config field, its default, and what changing it does | [Configuration reference](configuration.md) |
| Full output directory and column reference | [Output reference](outputs.md) |
| A worked example on public ENCODE mouse data | [Public mouse smoke test](PUBLIC_MOUSE.md) |
| Errors and troubleshooting | [FAQ](faq.md) |
