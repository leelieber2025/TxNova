# Tutorials

HTML on Read the Docs is pre-executed (tables and figures already there).

| If you want… | Open |
|--------------|------|
| FASTQ → STAR/HISAT2 BAM, strandedness, preflight | {doc}`t_prepare_bams` |
| Recover unannotated residual loci (main task) | {doc}`t_residual_catalog` |
| Rank gated loci, exon maps, fold viewer (GSE221720) | {doc}`t_gse221720_rank_and_fold` |
| 3Dmol fold viewer (example ORF) | {doc}`fold_viewer` |

**If you are new:** start with {doc}`t_prepare_bams` if you still have FASTQ,
then {doc}`../quickstart`, then {doc}`t_residual_catalog`. Output files:
{doc}`../outputs`.

## Run locally

```bash
pip install txnova pandas matplotlib
jupyter lab docs/tutorials/
```

| File | Used by |
|------|---------|
| `tests/fixtures/` | BAM preflight |
| `docs/tutorials/data/residual.tsv` | residual catalog (main-task notebook) |
| `docs/tutorials/data/gene_rank.tsv` | GSE221720 rank + fold |
| `docs/tutorials/data/candidates.tsv` | contrast-screen finals (bundled example) |
| `docs/tutorials/data/candidates.unnamed.tsv` | both-group (TPM) |

See `docs/tutorials/data/SOURCE.txt`. Residual IDs are one run's, not genes.

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item-card} Preparing input BAMs
:link: t_prepare_bams
:link-type: doc

STAR / HISAT2 recipes, strandedness, preflight on fixtures.
+++
Commands + tiny fixtures
:::

:::{grid-item-card} Residual catalog
:link: t_residual_catalog
:link-type: doc

Main task: harvest unannotated splices. Contrast is optional.
+++
Bundled residual.tsv · no BAM
:::

:::{grid-item-card} Rank and fold — GSE221720
:link: t_gse221720_rank_and_fold
:link-type: doc

Which gated residuals look gene-like, exon maps, TxNova 3Dmol viewer.
+++
Bundled TSVs · no BAM
:::

:::{grid-item-card} Fold viewer
:link: fold_viewer
:link-type: doc

Inline 3Dmol viewer TxNova writes for ORFs (pLDDT coloring).
+++
Example model
:::
::::

```{toctree}
:hidden: true
:maxdepth: 1

t_prepare_bams
t_residual_catalog
t_gse221720_rank_and_fold
fold_viewer
```
