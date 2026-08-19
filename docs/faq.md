# FAQ / Troubleshooting

First run: [Preparing input BAMs](tutorials/t_prepare_bams.ipynb) if you
still have FASTQ, then [Quickstart](quickstart.md). Main task:
[Residual catalog](tutorials/t_residual_catalog.ipynb). Every config field:
[Configuration reference](configuration.md). Output columns:
[Output reference](outputs.md).

## `txnova._core is not built`

This means you're on a source install and the Rust extension module didn't
compile — reinstall from your checkout:

```bash
pip install -e ".[dev]"
```

Check the install output for a Rust/`cargo` error — you need the Rust
toolchain installed (`rustup`) to build from source. If you installed with
plain `pip install txnova` (prebuilt wheel) and still see this, something
is wrong with the wheel for your platform — file an issue. See
[Installation](installation.md).

## Do I need extra tools on `PATH`?

`pip install txnova` is a prebuilt wheel. Residual splices are built inside
the process. Rust is only needed to develop the engine (see
[Installation](installation.md#development-install)).

## Preflight fails with a contig/`@SQ` message

```text
sample ctrl_1: BAM contig chr1 is not in FASTA .fai (BAM e.g. [...], FASTA e.g. [...])
```

or

```text
sample ctrl_1: contig chr1 LN=195471971 but FASTA .fai has 195154279
```

The BAM's `@SQ` header doesn't match `genome.fasta`'s `.fai` — either a
naming convention mismatch (`chr1` vs `1`) or a different genome build /
patch version. Names and lengths must match **literally**. Re-align against
the same FASTA you're pointing `genome.fasta` at, or fix the config to point
at the FASTA you actually aligned against.

## Preflight fails with "samples mix aligner families" or "mix library layouts"

Every sample in one run must come from the same aligner (all STAR or all
HISAT2) and the same layout (all paired-end or all single-end). Split mixed
sample sets into separate runs, or re-align the outliers.

## `group must be control|treat` / `strandedness must be unstranded|fr|rf`

`group` is optional. If the column is present, values must be exactly
`control` or `treat` — no `wt`/`ko`/`case` aliases. Omit the column for a
catalog-only run. Strandedness must be `unstranded`, `fr`, or `rf`. See
[Data preparation § figure out
strandedness](tutorials/t_prepare_bams.ipynb#5-figure-out-strandedness) if you're not
sure which of `fr`/`rf` your library is. The main task does not need both
groups: [Residual catalog](tutorials/t_residual_catalog.ipynb).

## `INFO txnova.orchestrator: DE skipped: need ≥2 control and ≥2 treat`

Not an error — this is TxNova telling you it turned DE off for this run.
PyDESeq2 needs replicates, so `de.enabled: true` in `config.yaml` only
takes effect when the sample sheet has ≥2 control and ≥2 treat rows;
below that, TxNova disables DE for you and logs this line instead of
failing. Every gate-passing locus then goes straight to the coding stage
(or straight to `candidates.tsv` if `coding.enabled: false` too) without a
DE filter. Add replicates if you want the DE filter to run; there's
nothing to fix if you don't. See [`de`](configuration.md#de).

## `candidates.tsv` is empty

An empty `candidates.tsv` is not the same as an empty run. Look at
`candidates/residual.tsv` first — that is the harvest catalog.

Without control and treat, `candidates.tsv` is the structure-pass
catalog; empty means every residual failed a structure gate.

With a contrast, look next at `candidates.unnamed.tsv` (also in control
by TPM) and `candidates.shared.tsv` (splice in both groups). Those two
are both-group unannotated structure, not a failed screen. See
[Output reference § the three tables](outputs.md#the-three-tables).

A contrast screen on a 2-vs-2 or 3-vs-3 design is often empty (see
[Public mouse smoke test](PUBLIC_MOUSE.md)). To see which stage cut the
**screen** count:

1. Open `report/report.md` (or `.html`) and read the **Funnel** section —
   it shows counts at every stage (merged transcripts → class `u` →
   structure-pass → gates → DE → final).
2. `candidates/candidates.gates.tsv` is the pre-DE view. If it has rows and
   `candidates.tsv` doesn't, the DE *filter* dropped them. `padj` NA with
   LFC ≥ `min_log2fc` is `de_status=low_count` and **stays**. Empty after
   gates means an evaluated non-hit (`padj` too large or LFC too small) or
   an unfitted row (`pvalue` also NA). Check `quantify/de.tsv`.
3. If `candidates.gates.tsv` is already empty, the structural/abundance
   gates in [`filters`](configuration.md#filters) are the ones to loosen
   first — most commonly `control_max_tpm`, `treat_detect_tpm`,
   `treat_median_tpm`, or the `discontinuity_*` group.
4. `candidates/residual.tsv` and `candidates/leak.tsv` are the splice
   harvest: cohort-recurrent CIGAR `N` junctions not in the annotation
   (`cohort=silent`, `shared`, or `cohort` when the sheet has no controls).
   A silent locus with strong `support_sum` but nothing in `candidates.tsv`
   failed a later gate — look at length, splice, distance, or (if you have
   a contrast) control TPM. Both-group splices land in
   `candidates.shared.tsv` when they pass structure gates.

## Unknown field / extra keys error from `config.yaml`

`config.yaml` is validated strictly (`extra="forbid"`) — a typo'd key, or a
key from an older TxNova version, raises an error naming the file instead of
being silently ignored. Compare against a freshly generated
`txnova init` file or [Configuration reference](configuration.md).

## What's the difference between `genome.annotation` and `genome.naming_annotation`?

`genome.annotation` is the **known-gene universe** — class `u` and leak
nearest genes are defined against it. Use a comprehensive annotation
(GENCODE comprehensive, not a cellranger-thin GTF); a thin one makes real
genes look intergenic.

`genome.naming_annotation` is optional. It fills `named_gene_name` /
`named_overlap` for class-`u` loci that `annotation` left blank, and residual
harvest unions it with `annotation` for the 200 nt same-strand knife (intron
span). It does not change class `u`. The 1 kb structure gate uses the run
`annotation` and the clipped locus. If `annotation` is already comprehensive,
you usually do not need it.

## `... has no gene or transcript rows; residual harvest would skip every gene filter`

Residual harvest needs `gene`/`transcript` (or `exon`, to build gene bodies
from) rows in `genome.annotation` and `genome.naming_annotation` to filter
new junctions against known gene bodies. A GTF with none of those rows —
empty file, wrong feature-type filter upstream, or a header-only stub — now
raises this error instead of silently harvesting with zero gene filters
(which would misclassify normal intragenic splices as intergenic). Point
`genome.annotation` at a comprehensive GTF that actually has gene/transcript
rows.

## Are locus IDs stable across runs?

No. `locus_id` (`RSDL.*` or a GENCODE id) is stable only inside one run.
A new assembler version or a different sample set can renumber residual
loci. Match candidates across runs by coordinates
(`chrom`/`start`/`end`/`strand`), not by ID.

## `ImportError: pydeseq2` or similar

`pydeseq2` is a required base dependency (not an optional extra) — if you
see an import error for it, the install is incomplete. Reinstall:

```bash
pip install txnova
# or, on a source install:
pip install -e ".[dev]"
```

## `coding.fold` / `coding.orphan` steps are slow or warn about network failures

Both make outbound calls to public APIs (AlphaFold DB, ESMFold, UniProt,
UCSC, EBI HMMER) and are on by default. They're fail-soft — a network error
is recorded as a warning, not a fatal error — but on a slow connection or an
offline cluster node they add real wall-clock time. Only the top 30
gene-like loci are folded. Set `coding.fold: false`
and/or `coding.orphan: false` in `config.yaml` to skip them. See
[Configuration reference § coding](configuration.md#coding).

## What does `threads: 0` actually do?

`0` means "use up to the CPU ceiling," but the number of BAM workers that
actually run in parallel at any moment is further capped by available
memory (`MemAvailable`) at run time — it's not a fixed pool size. Set an
explicit integer only if you want a hard ceiling below what memory would
otherwise allow (e.g. to share a machine with other jobs).

## A run I edited the config for didn't pick up my change

Each pipeline stage is skipped if its fingerprinted inputs (files + the
relevant slice of config) haven't changed since the last successful run —
you'll see `stamp hit <stage>` in the logs. If you changed something that
the fingerprint doesn't cover, or you just want a clean rebuild:

```bash
txnova run -c config.yaml --force
```

`--force` clears `output_dir/stamps/` and reruns every stage.

## `txnova report` says "no preflight.json"

`txnova report` rebuilds Markdown/HTML from an **existing** run's outputs —
it doesn't run the pipeline. Run `txnova run -c config.yaml` at least once
first.

## Coverage/structure gates behave strangely (e.g. everything looks non-canonical, or the funnel collapses at the structure step)

The most common cause is a wrong `strandedness` value in `samples.tsv` —
TxNova doesn't infer it from the BAM. See [Data preparation § figure out
strandedness](tutorials/t_prepare_bams.ipynb#5-figure-out-strandedness) to verify it
empirically before re-running.

## What does a `coding_label = coding` call actually mean?

The hexamer log-likelihood cleared `coding.hexamer_coding_min` (default
`0.0`), using the packaged table for `species` or your `hexamer_table`.
That threshold is the sign of the hexamer LLR, not CPAT's logistic cutoff.
Published CPAT cutoffs do not apply. `require_orf: true` is a separate gate:
a complete ORF of at least `min_orf_aa`.

## Mouse or human?

Only **mouse** and **human** are supported. Other species fail closed.

`species: auto` (default) reads the GTF (`#!genome-build`, then `ENSG…` /
`ENSMUSG…`), then `genome.assembly`. Hexamer table and UCSC tracks follow.
Set `species: mouse` or `species: human` to force it; that value must match
the GTF.

Human: GRCh38 FASTA + comprehensive GENCODE GTF. Mouse: GRCm39 + GENCODE
M39 comprehensive. Do not use a cellranger-thin GTF as the run annotation.

---

Still stuck? Check `output_dir/run.json` for the exact versions and config
used, and `output_dir/preflight.json` for the full preflight report.
