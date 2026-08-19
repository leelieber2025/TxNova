from __future__ import annotations

from pathlib import Path

import pytest

from txnova.config import (
    TxNovaConfig,
    load_config,
    peek_gtf_species,
    species_from_assembly,
)
from txnova.errors import TxNovaError


def test_unknown_key_forbidden(tmp_path: Path) -> None:
    p = tmp_path / "c.yaml"
    p.write_text(
        """
species: mouse
output_dir: ./out
threads: 0
genome:
  fasta: /a.fa
  annotation: /a.gtf
samples: s.tsv
not_a_real_key: 1
""",
        encoding="utf-8",
    )
    with pytest.raises(TxNovaError, match="invalid config"):
        load_config(p)


def test_defaults(tmp_path: Path) -> None:
    sheet = tmp_path / "s.tsv"
    sheet.write_text("sample_id\tbam\tgroup\tstrandedness\n", encoding="utf-8")
    p = tmp_path / "c.yaml"
    p.write_text(
        """
genome:
  fasta: /a.fa
  annotation: /a.gtf
samples: s.tsv
""",
        encoding="utf-8",
    )
    cfg = load_config(p)
    assert cfg.threads == 0
    assert cfg.genome.naming_annotation is None
    assert cfg.filters.class_ == "u"
    assert cfg.de.min_log2fc == 0.5
    assert cfg.samples == sheet.resolve()
    dumped = cfg.model_dump_canonical()
    assert dumped["filters"]["class"] == "u"


def test_model_direct() -> None:
    cfg = TxNovaConfig.model_validate(
        {
            "genome": {"fasta": "/a.fa", "annotation": "/a.gtf"},
            "samples": "/s.tsv",
        }
    )
    assert cfg.species == "auto"
    assert cfg.quantify.library_layout == "auto"
    assert cfg.coding.hexamer_table is None
    assert cfg.coding.fold is True
    assert cfg.coding.orphan is True
    assert cfg.filters.reject_bridging_junction is True
    assert cfg.filters.treat_min_detected_replicates == 3
    assert cfg.filters.control_max_tpm == 0.5
    assert cfg.filters.max_rmsk_frac == 0.1
    from txnova.config import packaged_hexamer_table, resolve_hexamer_table

    assert packaged_hexamer_table().is_file()
    assert packaged_hexamer_table("human").name == "Human_Hexamer.tsv"
    assert packaged_hexamer_table("human").is_file()
    assert resolve_hexamer_table(cfg).name == "Mouse_Hexamer.tsv"
    human = TxNovaConfig.model_validate(
        {
            "species": "human",
            "genome": {"fasta": "/a.fa", "annotation": "/a.gtf", "assembly": "GRCh38"},
            "samples": "/s.tsv",
        }
    )
    assert human.species == "human"
    assert resolve_hexamer_table(human).name == "Human_Hexamer.tsv"


def test_infer_species_from_gtf_and_assembly(tmp_path: Path) -> None:
    assert species_from_assembly("GRCh38") == "human"
    assert species_from_assembly("mm39") == "mouse"
    human_gtf = tmp_path / "h.gtf"
    human_gtf.write_text(
        "#!genome-build GRCh38.p14\n"
        'chr1\tHAVANA\tgene\t1\t10\t.\t+\t.\tgene_id "ENSG00000139618";\n',
        encoding="utf-8",
    )
    mouse_gtf = tmp_path / "m.gtf"
    mouse_gtf.write_text(
        "#!genome-build GRCm39\n"
        'chr1\tHAVANA\tgene\t1\t10\t.\t+\t.\tgene_id "ENSMUSG00000000001";\n',
        encoding="utf-8",
    )
    assert peek_gtf_species(human_gtf) == "human"
    assert peek_gtf_species(mouse_gtf) == "mouse"
    sheet = tmp_path / "s.tsv"
    sheet.write_text("sample_id\tbam\tgroup\tstrandedness\n", encoding="utf-8")
    cfg_path = tmp_path / "c.yaml"
    cfg_path.write_text(
        f"genome:\n  fasta: /a.fa\n  annotation: {human_gtf}\nsamples: s.tsv\n",
        encoding="utf-8",
    )
    cfg = load_config(cfg_path)
    assert cfg.species == "human"
    assert cfg.genome.assembly == "GRCh38"


def test_species_mismatch_gtf_vs_yaml(tmp_path: Path) -> None:
    gtf = tmp_path / "h.gtf"
    gtf.write_text("#!genome-build GRCh38\n", encoding="utf-8")
    sheet = tmp_path / "s.tsv"
    sheet.write_text("sample_id\tbam\tgroup\tstrandedness\n", encoding="utf-8")
    p = tmp_path / "c.yaml"
    p.write_text(
        f"species: mouse\ngenome:\n  fasta: /a.fa\n  annotation: {gtf}\n"
        f"  assembly: GRCh38\nsamples: s.tsv\n",
        encoding="utf-8",
    )
    with pytest.raises(TxNovaError, match="does not match the annotation GTF"):
        load_config(p)


def test_auto_species_human_gtf_overrides_mouse_assembly(tmp_path: Path) -> None:
    gtf = tmp_path / "h.gtf"
    gtf.write_text("#!genome-build GRCh38\n", encoding="utf-8")
    sheet = tmp_path / "s.tsv"
    sheet.write_text("sample_id\tbam\tgroup\tstrandedness\n", encoding="utf-8")
    p = tmp_path / "c.yaml"
    p.write_text(
        f"species: auto\ngenome:\n  fasta: /a.fa\n  annotation: {gtf}\n"
        f"  assembly: GRCm39\nsamples: s.tsv\n",
        encoding="utf-8",
    )
    cfg = load_config(p)
    assert cfg.species == "human"
    assert cfg.genome.assembly == "GRCh38"


def test_species_must_be_mouse_or_human() -> None:
    with pytest.raises(Exception):
        TxNovaConfig.model_validate(
            {
                "species": "zebrafish",
                "genome": {"fasta": "/a.fa", "annotation": "/a.gtf"},
                "samples": "/s.tsv",
            }
        )


def test_output_dir_resolves_against_config(tmp_path: Path) -> None:
    sheet = tmp_path / "s.tsv"
    sheet.write_text("sample_id\tbam\tgroup\tstrandedness\n", encoding="utf-8")
    p = tmp_path / "c.yaml"
    p.write_text(
        "genome:\n  fasta: /a.fa\n  annotation: /a.gtf\nsamples: s.tsv\noutput_dir: ./out\n",
        encoding="utf-8",
    )
    cfg = load_config(p)
    assert cfg.output_dir == (tmp_path / "out").resolve()


def test_homology_key_rejected(tmp_path: Path) -> None:
    sheet = tmp_path / "s.tsv"
    sheet.write_text("sample_id\tbam\tgroup\tstrandedness\n", encoding="utf-8")
    p = tmp_path / "c.yaml"
    p.write_text(
        "genome:\n  fasta: /a.fa\n  annotation: /a.gtf\nsamples: s.tsv\n"
        "coding:\n  homology: false\n",
        encoding="utf-8",
    )
    with pytest.raises(TxNovaError, match="invalid config"):
        load_config(p)


def test_stringtie_config_rejected() -> None:
    with pytest.raises(Exception):
        TxNovaConfig.model_validate(
            {
                "genome": {"fasta": "/a.fa", "annotation": "/a.gtf"},
                "samples": "/s.tsv",
                "stringtie": {"binary": "stringtie"},
            }
        )


def test_min_exons_one_allowed() -> None:
    cfg = TxNovaConfig.model_validate(
        {
            "genome": {"fasta": "/a.fa", "annotation": "/a.gtf"},
            "samples": "/s.tsv",
        }
    )
    assert cfg.filters.min_exons == 1
    one = TxNovaConfig.model_validate(
        {
            "genome": {"fasta": "/a.fa", "annotation": "/a.gtf"},
            "samples": "/s.tsv",
            "filters": {"min_exons": 1},
        }
    )
    assert one.filters.min_exons == 1


def test_min_exons_zero_rejected() -> None:
    with pytest.raises(Exception):
        TxNovaConfig.model_validate(
            {
                "genome": {"fasta": "/a.fa", "annotation": "/a.gtf"},
                "samples": "/s.tsv",
                "filters": {"min_exons": 0},
            }
        )
