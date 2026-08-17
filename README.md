# TxNova

[![PyPI version](https://img.shields.io/pypi/v/txnova.svg)](https://pypi.org/project/txnova/)
[![PyPI downloads](https://img.shields.io/pepy/dt/txnova.svg)](https://pepy.tech/project/txnova)
[![Python versions](https://img.shields.io/pypi/pyversions/txnova.svg)](https://pypi.org/project/txnova/)
[![Documentation Status](https://readthedocs.org/projects/txnova/badge/?version=latest)](https://txnova.readthedocs.io/en/latest/?badge=latest)
[![CI](https://github.com/leelieber2025/TxNova/actions/workflows/ci.yml/badge.svg)](https://github.com/leelieber2025/TxNova/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21970482.svg)](https://doi.org/10.5281/zenodo.21970482)

**TxNova** finds experimental-group-specific novel intergenic transcripts
from bulk RNA-seq BAMs. Treat-recurrent residual splices become locus
models. Class, counts, junctions, and bridges are recomputed from the BAM
and the annotation.

| Step | What it does |
|------|----------------|
| Harvest | Treat-recurrent CIGAR `N` junctions missing from the annotation → residual loci |
| Universe | Annotation + residual models, counted together |
| Gates | Structure, detection, optional DE |
| Tables | Treat-specific finals; both-group structure-pass (TPM and shared splice) |

Docs: [Read the Docs](https://txnova.readthedocs.io/en/latest/).

## What you need

- Python 3.10+
- Coordinate-sorted, indexed BAMs from STAR or HISAT2 (control and treat)
- Genome FASTA + `.fai`, and a comprehensive gene GTF (mouse default: GENCODE M39 / GRCm39)
- A YAML config and a sample sheet (`txnova init` writes starters)

## Install

```bash
pip install txnova
```

Wheels cover Linux and macOS. On Windows use WSL2. Building from source
needs a Rust toolchain; see
[Installation](https://txnova.readthedocs.io/en/latest/installation/).

## First run

```bash
txnova init -c config.yaml --samples samples.tsv
# edit paths, strandedness (rf / fr / unstranded), and sample rows
txnova preflight -c config.yaml
txnova run -c config.yaml
```

What to expect under `output_dir`:

- `candidates/candidates.tsv` — experimental-group-specific finals
- `candidates/candidates.unnamed.tsv` — structure-pass, also in control by TPM
- `candidates/candidates.shared.tsv` — structure-pass, splice in both groups
- `report/report.html` — start here
- `quantify/` — full-universe counts, TPM, and DE

The three candidate tables share the structural gates. They answer different
questions (induction vs unannotated structure in both groups). They are not
a declaration of new genes.

Next: [Quickstart](https://txnova.readthedocs.io/en/latest/quickstart/) ·
[Data preparation](https://txnova.readthedocs.io/en/latest/tutorials/t_prepare_bams/) ·
[Output reference](https://txnova.readthedocs.io/en/latest/outputs/) ·
[FAQ](https://txnova.readthedocs.io/en/latest/faq/)

## Status

**0.1.x.** Pin `txnova==0.1.3` in Methods. [Changelog](CHANGELOG.md).

## Citation

For the software, cite the Zenodo DOI above. Pin the installed
version in Methods (this tree is `txnova==0.1.3`). See `CITATION.cff`.

> Li, Z. TxNova. *Zenodo*. doi:10.5281/zenodo.21970482

PyPI: [https://pypi.org/project/txnova/](https://pypi.org/project/txnova/).

## License

Software: [Apache License 2.0](LICENSE). The packaged mouse hexamer table
comes from CPAT (Wang et al. 2013); see
[docs/license](https://txnova.readthedocs.io/en/latest/license/).

## Author

**Zhao Li (李钊)**  
Email: [leelieber@gmail.com](mailto:leelieber@gmail.com)
