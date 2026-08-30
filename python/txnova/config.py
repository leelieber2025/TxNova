from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Literal

import yaml
from pydantic import BaseModel, ConfigDict, Field, field_validator

from txnova.errors import TxNovaError

_HUMAN_ENSG = re.compile(r'gene_id "ENSG\d{11}')
_MOUSE_ENSG = re.compile(r'gene_id "ENSMUSG')


class _Strict(BaseModel):
    model_config = ConfigDict(extra="forbid")


class GenomeConfig(_Strict):
    fasta: Path
    annotation: Path
    annotation_source: str = "GENCODE"
    annotation_version: str = "M39"
    assembly: str = "GRCm39"
    # Optional comprehensive GTF (e.g. GENCODE M39) used only to name
    # class-u loci that the run annotation omitted. Not used for class u.
    naming_annotation: Path | None = None
    rmsk_bed: Path | None = None


class QuantifyConfig(_Strict):
    min_mapq: int = 10
    require_unique_nh: bool = False
    library_layout: Literal["auto", "paired", "single"] = "auto"
    skip_duplicate: Literal["auto", "always", "never"] = "auto"


class FiltersConfig(_Strict):
    class_: Literal["u"] = Field(default="u", alias="class")
    min_exons: int = Field(default=1, ge=1)
    require_canonical_splice: bool = True
    max_noncanonical_junction_fraction: float = Field(default=0.0, ge=0.0, le=1.0)
    min_nearest_same_strand_bp: int = Field(default=1000, ge=0)
    require_coverage_discontinuity: bool = True
    discontinuity_window_bp: int = Field(default=50, ge=1)
    discontinuity_valley_bp: int = Field(default=200, ge=1)
    discontinuity_valley_max_mean: float = Field(default=1.0, ge=0.0)
    discontinuity_ratio: float = Field(default=0.1, ge=0.0)
    discontinuity_min_treat_samples: int | None = Field(default=1, ge=1)
    reject_bridging_junction: bool = True
    bridge_min_reads: int = Field(default=2, ge=1)
    transcript_min_nt: int = Field(default=100, ge=1)
    control_max_tpm: float = Field(default=0.5, ge=0.0)
    treat_detect_tpm: float = Field(default=0.1, ge=0.0)
    treat_min_detected_replicates: int = Field(default=2, ge=1)
    treat_median_tpm: float = Field(default=0.5, ge=0.0)
    max_rmsk_frac: float = Field(default=0.1, ge=0.0, le=1.0)

    model_config = ConfigDict(extra="forbid", populate_by_name=True)


class DeConfig(_Strict):
    enabled: bool = True
    padj: float = Field(default=0.05, gt=0.0, le=1.0)
    min_log2fc: float = 0.5


class CodingConfig(_Strict):
    enabled: bool = True
    min_orf_aa: int = Field(default=50, ge=1)
    require_orf: bool = False
    hexamer_coding_min: float = 0.0
    hexamer_noncoding_max: float = 0.0
    hexamer_table: Path | None = None
    fold: bool = True
    orphan: bool = True


def normalize_species(value: str) -> str:
    s = (value or "").strip().lower()
    if s in {"auto", ""}:
        return "auto"
    if s in {"mouse", "mus musculus", "mmusculus", "mmu"}:
        return "mouse"
    if s in {"human", "homo sapiens", "hsapiens", "hsa"}:
        return "human"
    raise ValueError("species must be auto, mouse, or human (only mouse and human are supported)")


def species_from_assembly(assembly: str) -> str | None:
    a = (assembly or "").strip()
    if a in {"GRCh38", "hg38", "GRCh37", "hg19"}:
        return "human"
    if a in {"GRCm39", "mm39", "GRCm38", "mm10"}:
        return "mouse"
    return None


def peek_gtf_species(path: Path) -> str | None:
    """Read GENCODE headers / gene_id prefixes. Missing file → None."""
    if not path.is_file():
        return None
    genes = 0
    with path.open(encoding="utf-8", errors="replace") as fh:
        for i, line in enumerate(fh):
            if i > 400:
                break
            if line.startswith("#!"):
                low = line.lower()
                if "grch38" in low or "hg38" in low or "grch37" in low:
                    return "human"
                if "grcm39" in low or "mm39" in low or "grcm38" in low or "mm10" in low:
                    return "mouse"
                continue
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attrs = parts[8]
            if _MOUSE_ENSG.search(attrs):
                return "mouse"
            if _HUMAN_ENSG.search(attrs):
                return "human"
            genes += 1
            if genes >= 8:
                break
    return None


def default_assembly(species: str) -> str:
    return "GRCh38" if species == "human" else "GRCm39"


def apply_inferred_species(cfg: TxNovaConfig) -> TxNovaConfig:
    """Explicit mouse/human must match the GTF; auto reads GTF then assembly."""
    gtf_sp = peek_gtf_species(cfg.genome.annotation)
    asm_sp = species_from_assembly(cfg.genome.assembly)
    if cfg.species in {"mouse", "human"}:
        if gtf_sp and gtf_sp != cfg.species:
            raise TxNovaError(
                f"species: {cfg.species} does not match the annotation GTF ({gtf_sp})"
            )
    else:
        cfg.species = gtf_sp or asm_sp or "mouse"
    if species_from_assembly(cfg.genome.assembly) != cfg.species:
        cfg.genome.assembly = default_assembly(cfg.species)
    return cfg


class TxNovaConfig(_Strict):
    species: str = "auto"
    output_dir: Path = Path("./txnova_out")
    threads: int = 0
    genome: GenomeConfig
    samples: Path
    quantify: QuantifyConfig = Field(default_factory=QuantifyConfig)
    filters: FiltersConfig = Field(default_factory=FiltersConfig)
    de: DeConfig = Field(default_factory=DeConfig)
    coding: CodingConfig = Field(default_factory=CodingConfig)

    @field_validator("species")
    @classmethod
    def species_known(cls, v: str) -> str:
        return normalize_species(v)

    @field_validator("threads")
    @classmethod
    def threads_nonneg(cls, v: int) -> int:
        if v < 0:
            raise ValueError("threads must be >= 0 (0 = CPU ceiling; live workers follow memory)")
        return v

    def model_dump_canonical(self) -> dict[str, Any]:
        data = self.model_dump(mode="json", by_alias=True)
        return data


def load_config(path: Path) -> TxNovaConfig:
    path = path.resolve()
    if not path.is_file():
        raise TxNovaError(f"config not found: {path}")
    raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict):
        raise TxNovaError(f"config {path} must be a YAML mapping")
    try:
        cfg = TxNovaConfig.model_validate(raw)
    except Exception as e:
        raise TxNovaError(f"invalid config {path}: {e}") from e
    if not cfg.samples.is_absolute():
        cfg.samples = (path.parent / cfg.samples).resolve()
    if not cfg.genome.fasta.is_absolute():
        cfg.genome.fasta = (path.parent / cfg.genome.fasta).resolve()
    if not cfg.genome.annotation.is_absolute():
        cfg.genome.annotation = (path.parent / cfg.genome.annotation).resolve()
    if cfg.genome.naming_annotation is not None and not cfg.genome.naming_annotation.is_absolute():
        cfg.genome.naming_annotation = (path.parent / cfg.genome.naming_annotation).resolve()
    if not cfg.output_dir.is_absolute():
        cfg.output_dir = (path.parent / cfg.output_dir).resolve()
    if cfg.coding.hexamer_table is not None and not cfg.coding.hexamer_table.is_absolute():
        cfg.coding.hexamer_table = (path.parent / cfg.coding.hexamer_table).resolve()
    if cfg.genome.rmsk_bed is not None and not cfg.genome.rmsk_bed.is_absolute():
        cfg.genome.rmsk_bed = (path.parent / cfg.genome.rmsk_bed).resolve()
    if cfg.genome.rmsk_bed is not None and not cfg.genome.rmsk_bed.is_file():
        raise TxNovaError(f"genome.rmsk_bed not found: {cfg.genome.rmsk_bed}")
    return apply_inferred_species(cfg)


def packaged_hexamer_table(species: str = "mouse") -> Path:
    sp = normalize_species(species)
    if sp == "auto":
        sp = "mouse"
    name = "Human_Hexamer.tsv" if sp == "human" else "Mouse_Hexamer.tsv"
    return Path(__file__).resolve().parent / "data" / name


def resolve_hexamer_table(cfg: TxNovaConfig) -> Path:
    p = cfg.coding.hexamer_table
    if p is None:
        p = packaged_hexamer_table(cfg.species if cfg.species != "auto" else "mouse")
    if not p.is_file():
        raise TxNovaError(f"hexamer table not found: {p}")
    return p
