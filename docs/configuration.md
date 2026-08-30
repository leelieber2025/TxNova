# Configuration reference

`config.yaml` is validated strictly: every field is checked (types, ranges),
and **unknown fields raise an error** rather than being silently ignored.
`txnova init` writes a complete file with every default already filled in —
this page explains what each one means and when you'd want to change it.

Relative paths (`samples`, `genome.fasta`, `genome.annotation`,
`genome.naming_annotation`, `genome.rmsk_bed`, `coding.hexamer_table`) resolve
relative to the directory the config file is in, not your current working
directory. A configured `genome.rmsk_bed` that is missing is a hard error.

## Top level

| Field | Type | Default | Meaning |
|---|---|---|---|
| `species` | `auto`, `mouse`, or `human` | `auto` | Packaged hexamer table. `auto` reads the GTF then `genome.assembly`. `mouse` / `human` force that species (must match the GTF). Only mouse and human are supported. |
| `output_dir` | path | `./txnova_out` | Everything TxNova writes goes under here. See [Output reference](outputs.md). |
| `threads` | int ≥ 0 | `0` | `0` means "use up to the CPU ceiling"; the *live* number of parallel BAM workers is further capped by available memory at run time, not just this number. An explicit `N` is a hard ceiling, not a target. |
| `genome` | object | — | Required. See [`genome`](#genome). |
| `samples` | path | — | Required. Path to the sample sheet TSV. See [Quickstart § sample sheet](quickstart.md#3-fill-in-the-sample-sheet). |
| `quantify` | object | see below | See [`quantify`](#quantify). |
| `filters` | object | see below | See [`filters`](#filters). |
| `de` | object | see below | See [`de`](#de). |
| `coding` | object | see below | See [`coding`](#coding). |

## `genome`

| Field | Type | Default | Meaning |
|---|---|---|---|
| `fasta` | path | — | Required. Genome FASTA; must have a `.fai` index next to it. |
| `annotation` | path | — | Required. **This is the universe TxNova assembles against and classifies class `u` from.** Use a comprehensive annotation (e.g. GENCODE "comprehensive"), not a trimmed reference-package GTF — a thin annotation makes known genes look intergenic. |
| `annotation_source` | string | `GENCODE` | Provenance metadata, printed in the report. |
| `annotation_version` | string | `M39` | Provenance metadata, printed in the report. |
| `assembly` | string | `GRCm39` | UCSC conservation (`GRCm39`/`mm39` or `GRCh38`/`hg38`). If omitted and the GTF is human, set to `GRCh38`. |
| `naming_annotation` | path or `null` | `null` | Optional **second** GTF. It names class-`u` loci that `annotation` left blank (`named_gene_name` / `named_overlap`), and residual harvest unions it with `annotation` for the 200 nt same-strand knife (intron span). It does not change class `u`. The 1 kb structure gate uses the run `annotation`. |
| `rmsk_bed` | path or `null` | `null` | RepeatMasker BED (chrom start end name family). When set, `filters.max_rmsk_frac` is applied. Contig names are matched via `chrom_key` (`chr1` ≡ `1`). |

## `quantify`

Controls how reads are counted against the full universe (annotation +
residual loci), not just the final candidates.

| Field | Type | Default | Meaning |
|---|---|---|---|
| `min_mapq` | int | `10` | Minimum mapping quality (Phred; 10 ≈ 90%). Common counting floor. |
| `require_unique_nh` | bool | `false` | If `true`, drop reads unless `NH:i:1`. Multi-mappers are never down-weighted. |
| `library_layout` | `auto` \| `paired` \| `single` | `auto` | `auto` detects paired vs. single-end per BAM during preflight. Set explicitly only if you want preflight to *reject* a BAM that doesn't match. |
| `skip_duplicate` | `auto` \| `always` \| `never` | `auto` | Whether to skip reads flagged `0x400` (PCR/optical duplicate). `auto` counts all reads if preflight didn't see any duplicate-flagged reads in the BAM (i.e. you likely didn't mark duplicates), and skips them otherwise. |

## `filters`

Structural, detection, and abundance gates applied to class-`u`
(intergenic — no exon overlap, no gene body on either strand) merged loci.
A locus has to clear **all** of these to reach `candidates.gates.tsv` (before
DE and coding, which are separate stages).

| Field | Type | Default | Meaning |
|---|---|---|---|
| `class` | `u` (fixed) | `u` | Only `u` (fully intergenic) is supported. |
| `min_exons` | int ≥ 1 | `1` | Minimum exon count. `1` allows intronless loci through; the splice-canonicity gate below only applies when a locus actually has introns. |
| `require_canonical_splice` | bool | `true` | Reject loci whose scored introns are not **GT-AG, GC-AG, or AT-AC** (transcript strand). GT-AG/GC-AG are U2-type; AT-AC is the U12-type pair. |
| `max_noncanonical_junction_fraction` | float [0, 1] | `0.0` | Fraction of a locus's junctions allowed to be non-canonical before it's rejected (with `require_canonical_splice: true`). `0.0` = zero tolerance. |
| `min_nearest_same_strand_bp` | int ≥ 0 | `1000` | Minimum distance (bp) to the nearest same-strand annotated gene, measured on the **clipped residual locus**. Below this, a locus is treated as too close to be confidently intergenic. The default follows the ≥1 kb operational cutoff used for lincRNA catalogs, not typical gene–gene spacing (tens of kb). The 200 nt harvest knife is a separate cut and uses the intron span. |
| `require_coverage_discontinuity` | bool | `true` | Require a real coverage gap (valley) between the locus and its neighbor, not just annotation-based distance — guards against a locus that's actually the UTR/readthrough of an adjacent gene. |
| `discontinuity_window_bp` | int ≥ 1 | `50` | Window size (bp) used to sample coverage on either side of a candidate discontinuity. |
| `discontinuity_valley_bp` | int ≥ 1 | `200` | Length (bp) of the low-coverage valley required between locus and neighbor. |
| `discontinuity_valley_max_mean` | float ≥ 0 | `1.0` | Maximum mean coverage allowed inside the valley window for it to still count as a discontinuity. |
| `discontinuity_ratio` | float ≥ 0 | `0.1` | Valley coverage must drop to at most this fraction of flanking coverage. |
| `discontinuity_min_treat_samples` | int ≥ 1 or `null` | `1` | Minimum number of treat samples that must independently show the discontinuity. |
| `reject_bridging_junction` | bool | `true` | Reject a locus if there's a spliced read bridging it to the nearest same-strand gene — that's evidence it's actually part of that gene's transcript, not a separate locus. |
| `bridge_min_reads` | int ≥ 1 | `2` | Minimum spliced-bridge read count needed to trigger the rejection above. |
| `transcript_min_nt` | int ≥ 1 | `100` | Minimum spliced length (nt). Not the 200 nt lncRNA biotype cutoff — a shorter two-exon residual is still shown. |
| `control_max_tpm` | float ≥ 0 | `0.5` | Contrast only. Control **maximum** TPM must be **below** this for the **final** table. At or above → `candidates.unnamed.tsv`. Ignored when the sheet has no control+treat. |
| `max_rmsk_frac` | float 0–1 | `0.1` | Drop a locus if RepeatMasker covers this fraction of spliced length. Applied only when `genome.rmsk_bed` is set. |
| `treat_detect_tpm` | float ≥ 0 | `0.1` | Sample “detecting” the locus. Near the RNA-seq active/background floor. With a contrast, only treat samples count. |
| `treat_min_detected_replicates` | int ≥ 1 | `2` | Contrast only. Treat samples that must clear `treat_detect_tpm` for a locus to enter the final table. Default 2 so a first 2-vs-2 run can fill `candidates.tsv`. Harvest recurrence stays $r_{min}=2$ on the whole cohort. Preflight fails if this exceeds `n_treat`. |
| `treat_median_tpm` | float ≥ 0 | `0.5` | Contrast only. Treat median TPM. Half of the TPM 1 “on” line. |

If your final table is empty, `control_max_tpm` / `treat_detect_tpm` /
`treat_median_tpm` / `discontinuity_*` are the first knobs to loosen — see
the [FAQ](faq.md#candidatestsv-is-empty). Without control+treat, the final
table is the structure-pass residual set (no contrast filter).

## `de`

[PyDESeq2](https://pydeseq2.readthedocs.io/) fits the **full** locus-count
matrix (every universe locus) so size factors and FDR use the whole
transcriptome. The DE *filter* then keeps only gate-passing loci.

| Field | Type | Default | Meaning |
|---|---|---|---|
| `enabled` | bool | `true` | Contrast only. Skipped automatically unless the sheet has ≥2 control and ≥2 treat. Set `false` to skip DE even then. |
| `padj` | float (0, 1] | `0.05` | BH FDR cutoff (standard). |
| `min_log2fc` | float | `0.5` | Minimum **signed** log2 fold change (treat minus control; ≥1.4-fold). Only treat-up loci pass. |

A locus with `de_status = wald` met both `padj` and `min_log2fc`. A locus
with `de_status = low_count` has `padj` NA (DESeq2 independent filtering) but
a Wald `pvalue` and LFC ≥ `min_log2fc` — these **stay in the final table**.
Evaluated non-hits (padj too large, or LFC too small) are dropped. See
[Output reference](outputs.md#candidatestsv).

## `coding`

Optional ORF, coding-potential, structure, and evidence annotation on
gate-(and DE-)passing loci. Runs after `filters` and `de`.

| Field | Type | Default | Meaning |
|---|---|---|---|
| `enabled` | bool | `true` | Set `false` to skip the entire coding stage — `candidates.tsv` will still have structural/DE columns, just no ORF/coding/structure/conservation columns. |
| `min_orf_aa` | int ≥ 1 | `50` | Minimum ORF length (aa) to report. Conventional small-ORF boundary; CPAT searches from 75 nt. |
| `require_orf` | bool | `false` | If `true`, drop loci without a complete ORF of at least `min_orf_aa`. Off by default — many real intergenic loci are noncoding. |
| `hexamer_coding_min` | float | `0.0` | Hexamer log-likelihood above which a locus is called `coding_label = coding`. Published CPAT cutoffs do not apply. |
| `hexamer_noncoding_max` | float | `0.0` | Score at or below which a locus is called `noncoding`. Between `hexamer_noncoding_max` and `hexamer_coding_min` is `ambiguous`. No ORF at `min_orf_aa` is `no_orf`. |
| `hexamer_table` | path or `null` | `null` | `null` uses the packaged CPAT table for the resolved species (`Mouse_Hexamer.tsv` or `Human_Hexamer.tsv`). |
| `fold` | bool | `true` | Build 3D models for the **top 30** gene-like loci — AlphaFold DB for named loci, ESMFold for unnamed ones. **Requires internet access.** Network failures are recorded as warnings, not fatal. |
| `orphan` | bool | `true` | Same top 30: UCSC conservation (phyloP/phastCons) and EBI HMMER/Pfam. **Requires internet access.** Same fail-soft behavior. |

See [Installation § network access](installation.md#optional-network-access-for-structure-and-conservation)
if you need to disable `fold`/`orphan` for an offline run.

## Sample sheet (`samples.tsv`)

Not part of `config.yaml`, but referenced by its `samples:` field — see
[Quickstart § sample sheet](quickstart.md#3-fill-in-the-sample-sheet) for the
column reference.
