# Configuration reference

`config.yaml` is validated strictly: every field is checked (types, ranges),
and **unknown fields raise an error** rather than being silently ignored.
`txnova init` writes a complete file with every default already filled in —
this page explains what each one means and when you'd want to change it.

Relative paths (`samples`, `genome.fasta`, `genome.annotation`,
`genome.naming_annotation`, `coding.hexamer_table`) resolve relative to the
directory the config file is in, not your current working directory.

## Top level

| Field | Type | Default | Meaning |
|---|---|---|---|
| `species` | string | `mouse` | Provenance only. The packaged hexamer table is always mouse; this field does not switch tables. |
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
| `assembly` | string | `GRCm39` | Provenance metadata, printed in the report and used to tag structure/orphan lookups (e.g. UCSC track selection). |
| `naming_annotation` | path or `null` | `null` | Optional **second** GTF used only to *name* class-`u` loci that `annotation` doesn't have a gene body for. It is never used for class assignment or the distance gate — only `named_gene_name` / `named_overlap`. Residual harvest also uses it for the 200 nt same-strand knife. |

## `quantify`

Controls how reads are counted against the full universe (annotation +
residual loci), not just the final candidates.

| Field | Type | Default | Meaning |
|---|---|---|---|
| `min_mapq` | int | `10` | Minimum mapping quality for a read to be counted. |
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
| `require_canonical_splice` | bool | `true` | Reject loci whose scored introns are not **GT-AG or GC-AG** (transcript strand). AT-AC is not canonical here. |
| `max_noncanonical_junction_fraction` | float [0, 1] | `0.0` | Fraction of a locus's junctions allowed to be non-canonical before it's rejected (with `require_canonical_splice: true`). `0.0` = zero tolerance. |
| `min_nearest_same_strand_bp` | int ≥ 0 | `1000` | Minimum distance (bp) to the nearest same-strand annotated gene. Below this, a locus is treated as too close to be confidently intergenic. |
| `require_coverage_discontinuity` | bool | `true` | Require a real coverage gap (valley) between the locus and its neighbor, not just annotation-based distance — guards against a locus that's actually the UTR/readthrough of an adjacent gene. |
| `discontinuity_window_bp` | int ≥ 1 | `50` | Window size (bp) used to sample coverage on either side of a candidate discontinuity. |
| `discontinuity_valley_bp` | int ≥ 1 | `200` | Length (bp) of the low-coverage valley required between locus and neighbor. |
| `discontinuity_valley_max_mean` | float ≥ 0 | `1.0` | Maximum mean coverage allowed inside the valley window for it to still count as a discontinuity. |
| `discontinuity_ratio` | float ≥ 0 | `0.1` | Valley coverage must drop to at most this fraction of flanking coverage. |
| `discontinuity_min_treat_samples` | int ≥ 1 or `null` | `1` | Minimum number of treat samples that must independently show the discontinuity. |
| `reject_bridging_junction` | bool | `true` | Reject a locus if there's a spliced read bridging it to the nearest same-strand gene — that's evidence it's actually part of that gene's transcript, not a separate locus. |
| `bridge_min_reads` | int ≥ 1 | `2` | Minimum spliced-bridge read count needed to trigger the rejection above. |
| `transcript_min_nt` | int ≥ 1 | `100` | Minimum transcript length (nt) for a candidate. |
| `control_max_tpm` | float ≥ 0 | `1.0` | Control-group **maximum** TPM must be **below** this. At or above → `candidates.unnamed.tsv` (also in control), not the treat-specific table. |
| `treat_detect_tpm` | float ≥ 0 | `0.1` | TPM threshold above which a treat sample counts as "detecting" the locus. |
| `treat_min_detected_replicates` | int ≥ 1 | `2` | Minimum number of treat samples that must clear `treat_detect_tpm` ("recurrent in treat", not a one-off). Preflight also enforces this as a minimum sample count. |
| `treat_median_tpm` | float ≥ 0 | `0.5` | Median TPM across treat samples must be at least this. |

If your final table is empty, `control_max_tpm` / `treat_detect_tpm` /
`treat_median_tpm` / `discontinuity_*` are the first knobs to loosen — see
the [FAQ](faq.md#candidatestsv-is-empty).

## `de`

[PyDESeq2](https://pydeseq2.readthedocs.io/) fits the **full** locus-count
matrix (every universe locus) so size factors and FDR use the whole
transcriptome. The DE *filter* then keeps only gate-passing loci.

| Field | Type | Default | Meaning |
|---|---|---|---|
| `enabled` | bool | `true` | Set `false` to skip DE and take every gate-passing locus straight to coding/candidates (needs `de.enabled: false` if you have only 1 replicate per group, since DE needs ≥2 vs. ≥2 — preflight enforces this). |
| `padj` | float (0, 1] | `0.05` | BH-adjusted p-value cutoff. |
| `min_log2fc` | float | `0.5` | Minimum **signed** log2 fold change (treat minus control). Only treat-up loci pass. |

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
| `min_orf_aa` | int ≥ 1 | `50` | Minimum ORF length (amino acids) to be reported as `longest_orf_aa`; also the cutoff used by `require_orf`. |
| `require_orf` | bool | `false` | If `true`, drop loci without a complete ORF of at least `min_orf_aa`. Off by default — many real intergenic loci are noncoding. |
| `hexamer_coding_min` | float | `0.0` | Hexamer log-likelihood above which a locus is called `coding_label = coding`. Published CPAT cutoffs do not apply. |
| `hexamer_noncoding_max` | float | `0.0` | Score at or below which a locus is called `noncoding`. Between `hexamer_noncoding_max` and `hexamer_coding_min` is ambiguous. |
| `hexamer_table` | path or `null` | `null` | `null` uses the packaged mouse hexamer table (`python/txnova/data/Mouse_Hexamer.tsv`). Pass your own TSV to score a different species. |
| `fold` | bool | `true` | Build 3D structure models for predicted ORFs — AlphaFold DB for named loci, ESMFold for unnamed ones. **Requires internet access.** Network failures are recorded as warnings, not fatal. |
| `orphan` | bool | `true` | For loci with `named_overlap = none`, fetch UCSC conservation (phyloP/phastCons) and EBI HMMER/Pfam domain hits. **Requires internet access.** Same fail-soft behavior. |

See [Installation § network access](installation.md#optional-network-access-for-structure-and-conservation)
if you need to disable `fold`/`orphan` for an offline run.

## Sample sheet (`samples.tsv`)

Not part of `config.yaml`, but referenced by its `samples:` field — see
[Quickstart § sample sheet](quickstart.md#3-fill-in-the-sample-sheet) for the
column reference.
