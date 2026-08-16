# TxNova

[![Documentation](https://readthedocs.org/projects/txnova/badge/?version=latest)](https://txnova.readthedocs.io/)

Experimental-group-specific novel intergenic transcripts from bulk RNA-seq **BAMs**.

Treat-recurrent residual splices become locus models (`RSDL.*`). Class, counts,
junctions, and bridges are recomputed from the BAM and the annotation. The
output is a candidate table.

```bash
pip install txnova
txnova init -c config.yaml --samples samples.tsv
txnova preflight -c config.yaml
txnova run -c config.yaml
```

Wheels cover Linux and macOS. On Windows use WSL2. Input is a
coordinate-sorted, indexed BAM from `STAR` or `HISAT2`.

Outputs under `output_dir`:

- `candidates/candidates.tsv` — final candidate loci (after structure gates, optional DE, optional coding)
- `report/report.md` and `report/report.html`
- `quantify/` — full-universe counts / TPM / DE

Default coding table is the packaged mouse hexamer frequencies
(`coding.hexamer_table` overrides). The score is a hexamer log-likelihood.

A 2-vs-2 public mouse smoke (ENCODE BMDM ± Lipid A) is in
[docs/PUBLIC_MOUSE.md](docs/PUBLIC_MOUSE.md).

## Documentation

| Guide | Covers |
|---|---|
| [Docs home](docs/index.md) | Overview and pipeline diagram |
| [Installation](docs/installation.md) | `pip install`; source build only if you develop |
| [Quickstart](docs/quickstart.md) | First run, end to end |
| [Data preparation](docs/data_preparation.md) | FASTQ → TxNova-ready BAM, picking `strandedness` |
| [Configuration reference](docs/configuration.md) | Every `config.yaml` field and its default |
| [Output reference](docs/outputs.md) | Output directory layout and column meanings |
| [FAQ / Troubleshooting](docs/faq.md) | Preflight errors, empty candidate tables, common gotchas |
