from __future__ import annotations

from pathlib import Path

import pytest

from txnova.config import TxNovaConfig, load_config
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
    assert cfg.quantify.library_layout == "auto"
    assert cfg.coding.hexamer_table is None
    assert cfg.coding.fold is True
    assert cfg.coding.orphan is True
    assert cfg.filters.reject_bridging_junction is True
    from txnova.config import packaged_hexamer_table

    assert packaged_hexamer_table().is_file()


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
