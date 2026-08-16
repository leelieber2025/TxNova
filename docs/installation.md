# Installation

## Requirements

- **Python 3.10+**
- A coordinate-sorted, indexed BAM (STAR or HISAT2). TxNova does not align.

`pip install txnova` installs a **prebuilt wheel**. You do **not** need Rust,
Cargo, StringTie, samtools, or any other bioinformatics binary on `PATH`.

## Install

```bash
pip install txnova
```

Then:

```bash
txnova --version
python -c "import txnova._core as core; print(core.core_version())"
```

If `import txnova._core` raises `ImportError`, there is no wheel for your
platform yet — see [Development install](#development-install) (that path
needs a Rust toolchain) or file an issue.

Core Python dependencies come with the wheel: `typer`, `pydantic`, `pyyaml`,
`pandas`, `rich`, `jinja2`, `pydeseq2`, `numpy`, `scipy`.

## Optional: network access for structure and conservation

Two steps make outbound HTTPS calls. Both default **on**:

| Config flag | Calls | What it adds |
|---|---|---|
| `coding.fold: true` | AlphaFold DB, ESMFold, UniProt | 3D models for predicted ORFs |
| `coding.orphan: true` | UCSC, EBI HMMER/Pfam | Conservation and Pfam hits for unnamed loci |

Network failures are warnings and do **not** fail the run. On an offline
node, set both to `false` (see [Configuration](configuration.md#coding)).

## Next steps

| Step | Page |
|------|------|
| Run your first analysis | [Quickstart](quickstart.md) |
| Every config field | [Configuration reference](configuration.md) |
| What the pipeline writes | [Output reference](outputs.md) |

## Development install

Only for changing the code or running the test suite. This compiles the
Rust engine; you need `rustc` / `cargo` (`rustup`).

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
git clone https://github.com/leelieber2025/TxNova.git
cd TxNova
pip install -e ".[dev]"
pytest -q
```

`[dev]` adds `pytest`, `ruff`, `mypy`, and `maturin`. After editing Rust:
`maturin develop`.

Docs are on [txnova.readthedocs.io](https://txnova.readthedocs.io/).
Every push to `main` rebuilds the site. One GitHub secret is enough:

1. Read the Docs → account → [API tokens](https://readthedocs.org/accounts/tokens/) → create a token.
2. GitHub repo → Settings → Secrets → `RTD_TOKEN` = that token.

Do not click Import on the Read the Docs website. The `docs` workflow creates
the `txnova` project if it is missing and triggers a build.

Local preview: `pip install -e ".[docs]"` then `mkdocs serve`.

## Releasing wheels

Push a tag `v*` (for example `v0.1.0`). GitHub Actions builds manylinux and
macOS wheels plus an sdist and publishes them to PyPI. Users then get those
wheels from `pip install txnova`. Configure a PyPI **trusted publisher**
for this repo, workflow `release.yml`, environment `pypi` — no API token
in the repo.
