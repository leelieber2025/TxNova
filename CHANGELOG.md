# Changelog

Format: [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [0.1.5] - 2026-08-16

- Version bump (0.1.3 already used)

## [0.1.3] - 2026-08-16

### Added

- `candidates.shared.tsv` — structure-pass splices in both groups
- `gene_rank.tsv`; fold / Foldseek / orphan on top 30 only
- Tutorials (BAM prep, GSE221720 rank, 3Dmol viewer)
- `LICENSE` in the sdist

### Changed

- `min_nearest_same_strand_bp` default 500 (was 1000)
- Docs: three tables; Sphinx book theme

### Fixed

- 0.1.0 sdist rejected by PyPI (missing `LICENSE`)

## [0.1.0] - 2026-08-16

First release.

- CLI: `init`, `preflight`, `run`, `report`
- Residual harvest (treat-recurrent CIGAR `N` → RSDL loci)
- Class `u`; 200 nt same-strand knife
- Quantify, TPM, optional DE (`wald` / `low_count`)
- Structure gates; ORF + hexamer LLR + Fickett; optional fold
- Linux / macOS wheels. Windows: WSL2
