# Public mouse bulk smoke

End-to-end **contrast** smoke on public 2-vs-2 polyA RNA-seq (control and
treat both present). The main product is the residual catalog and does
not need both groups — see [Residual catalog](tutorials/t_residual_catalog.ipynb).
Large BAMs stay off git.

## Design

| | |
| --- | --- |
| Control | ENCODE [ENCSR614DLJ](https://www.encodeproject.org/experiments/ENCSR614DLJ/) untreated C57BL/6J BMDM, 2 reps |
| Treat | ENCODE [ENCSR614KOV](https://www.encodeproject.org/experiments/ENCSR614KOV/) Lipid A 100 ng/mL, 480 min, 2 reps |
| Library | polyA plus, reverse strand → TxNova `rf` |
| Align | STAR 2.7.2a, GRCm39 / GENCODE cellranger-2024-A |
| Assemble | TxNova residual assembler (annotation + leak splices) |

Input is BAM only. TxNova does not realign.

Other public mouse bulks (fracture; implant osteomyelitis) are tracked in the
GitHub repo, not on this site.

## Funnel magnitude (this run)

Counts, not IDs. Residual locus numbers (`RSDL.*`) can change across assembler
versions. The table is an **early smoke** (thinner annotation; distance gate
5 kb; unsigned `log2FC ≥ 2`). It is order of magnitude, not what a current
run prints.

Current software defaults: same-strand distance 1 kb; signed
`min_log2fc: 0.5` (treat-up only); `de_status=low_count` (`padj` NA, Wald
`pvalue` present, LFC ≥ 0.5) stays in the final table.

| step | n |
| --- | --- |
| merged transcripts | ~131 000 |
| merged loci | ~35 000 |
| class `u` transcripts | ~1 300 |
| all-`u` loci | ~1 200 |
| of which single-exon | ~88% |
| multi-exon + canonical splice + distance + coverage valley | ~80 |
| plus control-near-absent and treat-replicated | single digits |
| after pydeseq2 | same single digits on this 2-vs-2 |

Most loci in the universe are known genes (`overlap`). Treat-specific
intergenic loci are rare in Lipid A 8 h BMDM; a short table is expected.
Re-run with the current assembler before quoting IDs or exact counts.

## What to trust on a passing row

- `junction_support`: STAR CIGAR `N` on the representative’s own introns.
- `bridge_read_count`: spliced bridge to the nearest same-strand gene. Default gate rejects those.
- Full-universe TPM / DE.

A transcribed gap to a highly expressed same-strand neighbor can still pass if there is **no** spliced bridge. Inspect that class by hand (coverage in the gap, not just the 200 bp valley).

## Reproduce

1. Download the ENCODE FASTQs for the two experiments above.
2. Align with STAR to GRCm39 so `@SQ` contains the GTF seqnames.
3. Point `samples.tsv` at the four BAMs (`rf`) and `config.yaml` at the FASTA + GTF.
4. `txnova preflight -c config.yaml` then `txnova run -c config.yaml`.
5. Open `report/report.html`.

Do not commit BAM, FASTQ, or `genome.fa`.
