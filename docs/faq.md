# FAQ / Troubleshooting

First run: [Quickstart](quickstart.md). Every config field:
[Configuration reference](configuration.md). Output columns:
[Output reference](outputs.md). Preparing BAMs:
[Data preparation](tutorials/t_prepare_bams.ipynb).

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

These two columns accept exactly those literal strings — no `wt`/`ko`/`case`
aliases for `group`, and strandedness must be spelled `unstranded`, `fr`, or
`rf`. See [Data preparation § figure out
strandedness](tutorials/t_prepare_bams.ipynb#5-figure-out-strandedness) if you're not
sure which of `fr`/`rf` your library is.

## `de.enabled requires ≥2 control and ≥2 treat`

PyDESeq2 needs replicates. Either add samples, or set `de.enabled: false` in
`config.yaml` — every gate-passing locus then goes straight to the coding
stage (or straight to `candidates.tsv` if `coding.enabled: false` too)
without a DE filter. See [`de`](configuration.md#de).

## `candidates.tsv` is empty

An empty **final** table is not the same as an empty run. Look at
`candidates.unnamed.tsv` (structure-pass, also in control by TPM) and
`candidates.shared.tsv` (structure-pass, splice in both groups) first.
Those two are both-group unannotated structure, not induction. See
[Output reference § the three tables](outputs.md#the-three-tables).

Many 2-vs-2/3-vs-3 designs produce zero or a handful of intergenic,
treat-specific loci (see [Public mouse smoke test](PUBLIC_MOUSE.md)). To see
which stage cut the **final** count:

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
   harvest: treat-recurrent CIGAR `N` junctions not in the annotation
   (`cohort=silent` or `shared`). A silent locus with strong `treat_sum`
   but nothing in `candidates.tsv` failed a later gate — look at length,
   splice, distance, or control TPM. Both-group splices land in
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
harvest also uses it for the 200 nt same-strand knife. It does not change
class `u` or the 500 bp distance gate. If `annotation` is already comprehensive,
you usually do not need it.

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
`0.0`), using the packaged mouse table or your `hexamer_table`. Published
CPAT cutoffs do not apply. `require_orf: true` is a separate gate: a complete
ORF of at least `min_orf_aa`.

## Does TxNova support species other than mouse?

The packaged `coding.hexamer_table` (`python/txnova/data/Mouse_Hexamer.tsv`)
is mouse-specific. Everything upstream of coding (assembly, classification,
quantify, structure, gates, DE) is annotation-driven and species-agnostic —
point `genome.*` at a different species' FASTA/GTF and it works. For the
coding stage, either accept that the hexamer score won't be well-calibrated
for a non-mouse genome, set `coding.enabled: false`, or supply your own
`coding.hexamer_table` (see [Configuration reference §
coding](configuration.md#coding)).

---

Still stuck? Check `output_dir/run.json` for the exact versions and config
used, and `output_dir/preflight.json` for the full preflight report.
