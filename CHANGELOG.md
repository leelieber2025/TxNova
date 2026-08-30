# Changelog

Format: [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [0.1.10] - 2026-08-29

- Default `treat_min_detected_replicates` is 2 so a first 2-vs-2 contrast run
  can fill `candidates.tsv`. Preflight fails closed if the threshold exceeds
  `n_treat`, if the sheet mixes blank `group` with labeled groups, or if a
  configured `genome.rmsk_bed` is missing.
- Unstranded leak junctions take strand from the intron motif (GT-AG / GC-AG /
  AT-AC), then XS/ts; `+` and `-` are not clustered together. Canonical splice
  also scores AT-AC.
- `require_unique_nh` is honored in structure, leak, and quantify. Missing `NH`
  is unique (STAR/HISAT2), never treated as 0.
- RepeatMasker uses `chrom_key` (`chr1` ≡ `1`); rmsk is applied before the
  contrast split so high-rmsk control-expressed loci leave both unnamed and
  shared. `gates` stamps fingerprint the BED (including null); residual stamps
  fingerprint FASTA `.fai`.
- Quantify requires a matching CIGAR `N` for spliced loci. Leak harvest omits
  annotation introns (no `assembled_u` class-hardcode). `annotate_leak` joins
  harvest junctions to residual `intron_structure`.
- DESeq uses only labeled control/treat columns. Paired-end BAM EOF flushes
  remaining pending mates as dropped fragments (leak and quantify).
- `coding_label` is `coding` / `noncoding` / `ambiguous` / `no_orf`.
  `min_orf_aa` is the reporting floor. `chrom_class("1")` and `("Y")` are
  primary. Shared-harvest IDs are excluded from finals.
- Shared `USER_AGENT = txnova/{version}`; Foldseek uses urllib. One shared
  FASTA `.fai` parser; a bad contig skips that locus instead of aborting the
  run. Preflight always sample-checks BAM coordinate order.

## [0.1.9] - 2026-08-23

- Fixed: the Folding request User-Agent now reports the installed package
  version instead of a stale hard-coded value.
- Removed unused Rust fields and redundant conversions; normal Rust builds are
  warning-free.
- Added repository, homepage, and documentation metadata to the Rust crate.
- Updated Python package license metadata to the current SPDX form.
- Docs: corrected user-facing wording and standardized American English.

## [0.1.8] - 2026-08-22

- Fixed: unstranded (`strandedness: unstranded`) loci had their canonical-splice check and their ORF/coding-score scan silently read on the plus strand only, so a true minus-strand GT-AG/GC-AG intron scored as non-canonical and some loci were translated on the wrong strand. `splice::infer_unstranded_plus` now scores each intron's donor/acceptor motif in both directions and picks the better-supported strand once per transcript; both `splice_features` (the canonical-splice gate) and `orf::splice_seq` (ORF finding / coding score) use that call instead of each defaulting to plus independently. Only unstranded samples are affected.

## [0.1.7] - 2026-08-19

- `residual.tsv` `nearest_gene_name` / `nearest_distance_bp` are the distance of the residual locus (the interval the 1 kb structure gate uses), not of the nearest member intron. Stranded loci use same-strand distance; unstranded loci use either-strand. The column is rewritten after terminal clip and after the coverage walk.

## [0.1.6] - 2026-08-18

- Residual harvest is cohort-recurrent (any ≥2 samples), not treat-recurrent. Sample `group` is optional. Control-versus-treat TPM gates and DE run only when both groups are present; DE is skipped unless each side has ≥2 samples.
- Docs: main-task tutorial `docs/tutorials/t_residual_catalog.ipynb` (catalog first; contrast optional). User-facing pages no longer require control and treat.
- Default `treat_min_detected_replicates` is 3 (locus detection gate). Junction harvest stays at r_min = 2.
- Default `control_max_tpm` is 0.5. Optional `genome.rmsk_bed` + `filters.max_rmsk_frac` (0.1) drops repeat-heavy residual models.
- Residual harvest keeps a locus when the intron chain is intact (`n_exons ≥ 2`), even if a terminal stays a 30 nt stub after the coverage walk. Dropping the whole chain was discarding leak-supported genes (mask-100: Disp2, Kprp, Krt2, Slc12a8).
- Do not chain residual junctions across a remaining gene body. The old merge-then-drop-span killed both sides (mask-1000: Dpf2 / Cdc42ep2 across Gm42067).
- GTF with no `gene` rows: build gene bodies from transcript/exon; missing or empty annotation raises instead of silently skipping every harvest knife.
- `chrom_key` maps chrM/M/MT to `MT`.
- Residual terminals clip to FASTA `.fai` length (no exon past contig end).
- Quantify indexes exons by contig name per BAM (not the first sample's tid).
- FASTA fetch reads a block instead of one seek per base.
- Unknown class codes and junction_support length mismatches raise `TxNovaError`.
- Terminal walk window includes the last `max_terminal_nt` step. Duplicate flag counts only mapped primary reads.

## [0.1.5] - 2026-08-17

- Version bump (0.1.3 already used)
- `min_nearest_same_strand_bp` default 1000 (lincRNA ≥1 kb cutoff)
- `species: auto` (default) infers mouse/human from the GTF; `mouse` / `human` override. Only those two are supported.

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
