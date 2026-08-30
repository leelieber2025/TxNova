from __future__ import annotations

from pathlib import Path

from txnova.errors import TxNovaError

CONFIG_TEMPLATE = """\
species: auto                  # auto | mouse | human. Only mouse and human are supported.
output_dir: ./txnova_out
threads: 0                 # 0 = CPU ceiling; live BAM workers follow MemAvailable. Explicit N is a ceiling.

genome:
  fasta: /path/to/GRCm39.primary_assembly.genome.fa
  annotation: /path/to/gencode.vM39.primary_assembly.annotation.gtf
  annotation_source: GENCODE
  annotation_version: M39              # human: e.g. v47
  assembly: GRCm39                     # human: GRCh38
  naming_annotation: null          # optional comprehensive GTF to name class-u loci the run annotation omitted

samples: samples.tsv

quantify:
  min_mapq: 10
  require_unique_nh: false
  library_layout: auto
  skip_duplicate: auto

filters:
  class: u
  min_exons: 1              # intronless genes are allowed; splice gate only if n_introns>0
  require_canonical_splice: true
  max_noncanonical_junction_fraction: 0.0
  min_nearest_same_strand_bp: 1000
  require_coverage_discontinuity: true
  discontinuity_window_bp: 50
  discontinuity_valley_bp: 200
  discontinuity_valley_max_mean: 1.0
  discontinuity_ratio: 0.1
  discontinuity_min_treat_samples: 1
  reject_bridging_junction: true
  bridge_min_reads: 2
  transcript_min_nt: 100
  control_max_tpm: 0.5
  treat_detect_tpm: 0.1
  treat_min_detected_replicates: 2
  treat_median_tpm: 0.5
  max_rmsk_frac: 0.1

de:
  enabled: true
  padj: 0.05
  min_log2fc: 0.5

coding:
  enabled: true
  min_orf_aa: 50
  require_orf: false
  hexamer_coding_min: 0.0
  hexamer_noncoding_max: 0.0
  hexamer_table: null      # null = packaged table for species (mouse or human)
  fold: true               # AlphaFold DB if named; ESMFold if unnamed ORF
  orphan: true             # named_overlap=none: UCSC phyloP/phastCons + HMMER/Pfam
"""

SAMPLES_TEMPLATE = """\
sample_id	bam	group	strandedness	replicate
ctrl_1	/data/ctrl_1.bam	control	rf	1
ctrl_2	/data/ctrl_2.bam	control	rf	2
treat_1	/data/treat_1.bam	treat	rf	1
treat_2	/data/treat_2.bam	treat	rf	2
"""


def write_init_files(config_path: Path, samples_path: Path) -> None:
    if config_path.exists():
        raise TxNovaError(f"refusing to overwrite {config_path}")
    if samples_path.exists():
        raise TxNovaError(f"refusing to overwrite {samples_path}")
    config_path.parent.mkdir(parents=True, exist_ok=True)
    samples_path.parent.mkdir(parents=True, exist_ok=True)
    cfg_dir = config_path.parent.resolve()
    samples_resolved = samples_path.resolve()
    try:
        samples_ref = str(samples_resolved.relative_to(cfg_dir))
    except ValueError:
        samples_ref = str(samples_resolved)
    body = CONFIG_TEMPLATE.replace("samples: samples.tsv", f"samples: {samples_ref}")
    config_path.write_text(body, encoding="utf-8")
    samples_path.write_text(SAMPLES_TEMPLATE, encoding="utf-8")
