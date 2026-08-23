# TxNova Documentation

[![PyPI version](https://img.shields.io/pypi/v/txnova.svg)](https://pypi.org/project/txnova/)
[![PyPI downloads](https://img.shields.io/pepy/dt/txnova.svg)](https://pepy.tech/project/txnova)
[![CI](https://github.com/leelieber2025/TxNova/actions/workflows/ci.yml/badge.svg)](https://github.com/leelieber2025/TxNova/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](https://github.com/leelieber2025/TxNova/blob/main/LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21970482.svg)](https://doi.org/10.5281/zenodo.21970482)

## What TxNova does

**TxNova** — short for **Transcript Nova** — recovers **unannotated spliced
residual loci** from existing bulk RNA-seq BAMs. Cohort-recurrent residual
splices become locus models. Class, counts, junctions, and bridges are
recomputed from the BAM and the annotation. A control-versus-treat filter
is optional.

| Step | What it does |
|------|----------------|
| Harvest | Cohort-recurrent CIGAR `N` junctions missing from the annotation → residual loci |
| Universe | Annotation + residual models, counted together |
| Gates | Structure always; treat-detected / control-silent only if both groups are present |
| Tables | Residual catalog; optional contrast finals; both-group tables when a contrast exists |

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
    │ 2. Leak + residual     │  cohort-recurrent CIGAR N (silent / shared / cohort)
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
    │ 6. Coding              │  ORF, hexamer LLR, Fickett; fold top 30
    └──────────────────────┘
```

Per-sample GTFs are stubs. Leak and quantify scan the BAMs. The primary
output is the residual catalog. When the sheet has both control and
treat, the run also writes a treat-detected / control-silent screen plus
unnamed (control TPM) and shared (splice in both groups).

## Where to go

| Goal | Page |
|------|------|
| Install | {doc}`installation` |
| Turn FASTQ into a TxNova-ready BAM | {doc}`tutorials/t_prepare_bams` |
| First run after you have BAMs | {doc}`quickstart` |
| Residual catalog (main task) | {doc}`tutorials/t_residual_catalog` |
| Every `config.yaml` field | {doc}`configuration` |
| What each output file and column means | {doc}`outputs` |
| Empty tables and preflight errors | {doc}`faq` |
| Ranked loci, exon maps, fold viewer | {doc}`tutorials/index` |
| Public mouse smoke (ENCODE BMDM ± Lipid A) | {doc}`PUBLIC_MOUSE` |
| Changelog | {doc}`changelog` |

### A sensible path

1. Install: `pip install txnova`.
2. If you still have FASTQ, align first ({doc}`tutorials/t_prepare_bams`).
3. Follow {doc}`quickstart`.
4. Open `report/report.html`, then `candidates/residual.tsv` and the
   tables in {doc}`outputs`.
5. Main-task walkthrough: {doc}`tutorials/t_residual_catalog`.

### Default call

```bash
txnova init -c config.yaml --samples samples.tsv
# edit paths, strandedness (rf / fr / unstranded), and sample rows
txnova preflight -c config.yaml
txnova run -c config.yaml
```

`txnova preflight` is cheap and catches mismatched contigs, missing indexes,
mixed aligners, and bad `group` / `strandedness` values. `txnova run` runs
preflight again. `txnova report` rebuilds Markdown/HTML from an existing run.

You need coordinate-sorted, indexed BAMs from STAR or HISAT2 (not a mix),
at least two samples (control and treat are optional; DE needs ≥2 of each),
a FASTA with `.fai`, and a comprehensive GTF whose seqnames match the BAM
`@SQ` lines.
Mouse (GENCODE M39 / GRCm39) and human (GENCODE / GRCh38) are supported.
`species: auto` infers which one; you can set `mouse` or `human`.

:::{note}
TxNova is **0.1.x** (mouse and human). The documented interface is the CLI
(`init`, `preflight`, `run`, `report`) and `config.yaml`. Pin the installed
version in Methods (`txnova==0.1.9` for this tree). Cite
[doi:10.5281/zenodo.21970482](https://doi.org/10.5281/zenodo.21970482).
:::

::::{grid} 1 2 3 3
:gutter: 2

:::{grid-item-card} Installation {octicon}`plug;1em;`
:link: installation
:link-type: doc

`pip install txnova`; source build if you develop.
:::

:::{grid-item-card} Preparing input BAMs {octicon}`database;1em;`
:link: tutorials/t_prepare_bams
:link-type: doc

STAR / HISAT2 recipes; preflight on fixtures.
:::

:::{grid-item-card} Quickstart {octicon}`rocket;1em;`
:link: quickstart
:link-type: doc

First run after you have BAMs.
:::

:::{grid-item-card} Configuration {octicon}`gear;1em;`
:link: configuration
:link-type: doc

Every YAML field and default.
:::

:::{grid-item-card} Outputs {octicon}`table;1em;`
:link: outputs
:link-type: doc

The three tables and column meanings.
:::

:::{grid-item-card} FAQ {octicon}`question;1em;`
:link: faq
:link-type: doc

Preflight errors and empty tables.
:::

:::{grid-item-card} Tutorials {octicon}`play;1em;`
:link: tutorials/index
:link-type: doc

Prepare BAMs first; then the residual catalog.
:::

:::{grid-item-card} Public mouse {octicon}`beaker;1em;`
:link: PUBLIC_MOUSE
:link-type: doc

ENCODE BMDM smoke test.
:::

:::{grid-item-card} License {octicon}`law;1em;`
:link: license
:link-type: doc

Apache-2.0 and hexamer table.
:::

:::{grid-item-card} GitHub {octicon}`mark-github;1em;`
:link: https://github.com/leelieber2025/TxNova

Source and issues.
:::
::::

```{toctree}
:hidden: true
:maxdepth: 1
:titlesonly: true

installation
quickstart
tutorials/index
configuration
outputs
faq
PUBLIC_MOUSE
changelog
license
GitHub <https://github.com/leelieber2025/TxNova>
```
