# Output reference

Everything a run produces lives under `output_dir` (default `./txnova_out`,
set via `config.yaml`'s top-level `output_dir`).

```text
output_dir/
├── preflight.json               # preflight check results
├── run.json                     # provenance: versions, config, metrics, warnings
├── stamps/                      # per-stage fingerprints, used to skip unchanged work
├── assembly/
│   ├── <sample_id>.gtf          # stub (no per-sample reconstruct)
│   ├── merged.gtf               # gene/transcript/exon rows from the annotation
│   └── universe.gtf             # merged.gtf + residual.gtf — the quantify universe
├── classify/
│   ├── transcripts.class.tsv    # class code per universe transcript
│   ├── representatives.tsv      # one representative transcript per locus
│   └── locus.class.tsv          # locus-level class (all isoforms of a gene_id)
├── quantify/
│   ├── locus_tpm.tsv            # TPM, every locus in the universe
│   ├── locus_counts.tsv         # raw counts, every locus in the universe
│   └── de.tsv                   # full PyDESeq2 results table (if de.enabled)
├── structure/
│   └── structure.features.tsv   # splice canonicity, coverage, bridging, distances
├── candidates/
│   ├── leak.tsv                 # cohort-recurrent splices missing from the annotation (silent / shared / cohort)
│   ├── residual.tsv             # locus hypotheses built from leaked splices
│   ├── residual.gtf
│   ├── residual.stats.json      # n_loci / n_degenerate
│   ├── candidates.gates.tsv     # after structural gates, before DE
│   ├── candidates.unnamed.tsv   # structure-pass loci that are also seen in control (TPM)
│   ├── candidates.shared.tsv    # structure-pass loci whose splices are also in control
│   ├── candidates.de.tsv        # after DE (or a passthrough copy if de.enabled: false)
│   ├── candidates.reps.tsv      # representative transcript per candidate locus
│   ├── orfs.tsv                 # ORF scan output
│   ├── peptides.fa              # translated longest ORFs
│   ├── orphan.tsv               # conservation + Pfam for named_overlap=none loci
│   ├── gene_rank.tsv            # structure-pass loci ranked by how gene-like they look
│   ├── transcripts.tsv          # all merged transcripts belonging to final candidate loci
│   └── candidates.tsv           # ← the final table
├── fold/
│   ├── index.html               # 3D structure viewer (3Dmol.js)
│   ├── structures.tsv
│   └── function.tsv             # UniProt (named) / Foldseek fold matches (unnamed)
└── report/
    ├── report.md
    └── report.html              # ← start here
```

## Start here: `report/report.html`

The report is a self-contained walkthrough of the run: parameters used, the
sample sheet, preflight results, a **funnel table** (how many transcripts/loci
survived each stage), and the three structure-pass tables plus gene-rank,
orphan, leak, and residual. `report/report.md` has the same content as plain
Markdown. Regenerate either without rerunning the pipeline:

```bash
txnova report -c config.yaml
```

## The three tables

A run writes **three** structure-pass views. They share length, distance,
splice, valley, and bridge gates. They are not three snapshots of the same
filter.

| File | Who is in it | What it means |
|---|---|---|
| `candidates/residual.tsv` | Clustered harvest models | The residual catalog. Always written. |
| `candidates/candidates.tsv` | Structure pass. With both groups: also control-silent and treat-detected (and DE if it ran) | Catalogue without a contrast; treat-detected / control-silent **screen** with a contrast. |
| `candidates/candidates.unnamed.tsv` | Structure pass, control max TPM at or above the gate | Interval also transcribed in control. Contrast only. |
| `candidates/candidates.shared.tsv` | Structure pass, and a harvest junction is `cohort=shared` | The splice is used in both groups. Contrast only. Not a DE call. |

A locus can appear in unnamed and shared at once. It cannot appear in both
the final table and unnamed (those two are split by the control TPM gate).
`candidates.gates.tsv` / `candidates.de.tsv` are earlier snapshots of the
**final** path only.

## `candidates.tsv`

The final table — structure-pass residual loci. When the sheet has control and treat, this is the treat-detected / control-silent subset. Without a contrast it is the residual catalog.

| Column | Meaning |
|---|---|
| `locus_id` | Stable ID for this locus within the run (not a persistent cross-run gene ID — see the [FAQ](faq.md#are-locus-ids-stable-across-runs)). |
| `locus_coord`, `chrom`, `start`, `end`, `strand` | Genomic coordinates. |
| `n_exons`, `length_nt`, `exon_structure` | Structure of the representative transcript. |
| `class` | Always `u` in the final table (fully intergenic — no exon overlap, no gene body on either strand). |
| `nearest_gene_id` / `nearest_gene_name` / `nearest_distance_bp` / `nearest_strand` | Nearest **same-strand** annotated gene. |
| `nearest_any_gene_id` / `nearest_any_gene_name` / `nearest_any_distance_bp` / `nearest_any_strand` | Nearest annotated gene on **either** strand. |
| `named_gene_name` / `named_gene_id` / `named_gene_type` / `named_overlap` | Filled from `genome.naming_annotation` if configured; `named_overlap = none` means neither annotation could name it (these loci feed `orphan.tsv`). |
| `canonical_splice_fraction` | Fraction of this locus's junctions that are GT-AG or GC-AG. |
| `coverage_discontinuity` / `coverage_valley_mean` / `coverage_gap_mean` | Evidence for a real coverage gap between this locus and its neighbor. |
| `control_max_tpm` / `treat_median_tpm` / `treat_n_detected` | Abundance gate inputs — see [`filters`](configuration.md#filters). |
| `junction_support` / `junction_support_min` | BAM CIGAR `N`-based splice support. |
| `bridge_read_count` / `has_bridge` | Spliced reads bridging to the nearest same-strand gene — `reject_bridging_junction` gates on this. |
| `gates_passed` | Which structural/detection gates this locus cleared. |
| `representative_transcript_id` | Which merged-transcript isoform this row's structure/coding columns describe. |
| `assembled_in` | Unused with the residual assembler (per-sample GTFs are stubs). |
| `baseMean`, `log2FoldChange`, `pvalue`, `padj` | PyDESeq2 output (present when `de.enabled: true`). |
| `de_status` | `wald` = `padj` and `min_log2fc` both met. `low_count` = DESeq2 left `padj` as NA (independent filtering) but Wald `pvalue` exists and LFC meets `min_log2fc` — these **stay** in the final table. |
| `de_pass` | Boolean; whether the row is in the final table because of DE (mirrors `de_status`). |
| `longest_orf_aa`, `orf_complete` | Longest ORF found and whether it has both start and stop codons. |
| `coding_score`, `fickett_score` | Hexamer log-likelihood and Fickett score. Published CPAT cutoffs do not apply. |
| `coding_label` | `coding` / `noncoding` / ambiguous, from `coding.hexamer_coding_min` / `hexamer_noncoding_max` in [`coding`](configuration.md#coding). |

Per-sample columns (`<sample_id>_junction_support`, `<sample_id>_bridge_read_count`,
`<sample_id>_tpm`, `<sample_id>_count`) are also included so you can inspect
per-replicate support for any row by hand.

## `candidates/gene_rank.tsv`

Structure-pass loci from the **final table ∪ unnamed** (not shared-only
rows that failed those two), ranked by how gene-like they look — ORF /
hexamer / Fickett plus GENCODE-style checks and penalties from manual
review (unplaced contigs, zero-support "introns," retroviral gag/pol/MLV,
intronless copies of known proteins). Stage 1 scores every row offline.
Stage 2 fetches phyloP / Pfam on the top **40**, then re-sorts; the report
shows the top **30**. Use this when `candidates.tsv` is short; unnamed
loci are in the rank on purpose.

## `candidates/orphan.tsv`

Conservation (`phylop_mean`, `phylop_frac_pos`, `phastcons_mean` — UCSC
tracks on exons, not PhyloCSF) and Pfam domain hits (EBI HMMER `hmmscan`,
`pfam_name` / `pfam_evalue` / `pfam_desc`) for loci where neither annotation
could supply a gene name (`named_overlap = none`).

## `candidates/leak.tsv` and `candidates/residual.tsv`

**`leak.tsv`** lists cohort-recurrent spliced junctions from the BAMs
(CIGAR `N`) that are absent from the annotation (`status = unassembled`).
`cohort=silent` has no control support; `cohort=shared` is also in control;
`cohort=cohort` is used when the sheet has no control samples.
Junctions inside known gene bodies are omitted.

**`residual.tsv`** clusters those junctions into locus hypotheses (shared
splice site, or a 30–20 kb constitutive exon between adjacent introns),
clips terminals off gene bodies, and writes them into `assembly/universe.gtf`.
They then take the same quantify → structure → gates → DE path as annotated
genes. `cohort` on each locus is `shared` if any member junction was in control.

See [The three tables](#the-three-tables) for how `candidates.shared.tsv`
and `candidates.unnamed.tsv` differ from the final table.

## `fold/`

3D structure models for the **top 30** gene-like loci (`gene_rank.tsv`) —
AlphaFold DB lookups for named loci, ESMFold predictions for unnamed peptides.
Open `fold/index.html` for an interactive viewer (colored by pLDDT confidence);
`fold/function.tsv` adds UniProt curated function (named loci) or Foldseek
fold-similarity hits (unnamed loci). Only produced when `coding.enabled` and
`coding.fold` are both `true`. Other loci are not folded.

## `run.json`

Machine-readable provenance: assembler version (`txnova-residual-*`), Rust
core, full config, per-stage metrics, and warnings.

## `stamps/`

Per-stage fingerprints (inputs + config hash) that `txnova run` checks before
redoing work — if nothing relevant changed since the last run, that stage is
skipped and you'll see `stamp hit <stage>` in the logs. Delete this directory,
or pass `--force` to `txnova run`, to force a full rebuild.
