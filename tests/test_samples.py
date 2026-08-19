from __future__ import annotations

from pathlib import Path

import pytest

from txnova.errors import TxNovaError
from txnova.samples import load_samples


def test_rejects_alias_group(tmp_path: Path) -> None:
    p = tmp_path / "s.tsv"
    p.write_text(
        "sample_id\tbam\tgroup\tstrandedness\na\ta.bam\tcase\trf\n",
        encoding="utf-8",
    )
    with pytest.raises(TxNovaError, match="control\\|treat"):
        load_samples(p)


def test_rejects_extra_column(tmp_path: Path) -> None:
    p = tmp_path / "s.tsv"
    p.write_text(
        "sample_id\tbam\tgroup\tstrandedness\tnote\na\ta.bam\tcontrol\trf\tx\n",
        encoding="utf-8",
    )
    with pytest.raises(TxNovaError, match="unknown columns"):
        load_samples(p)


def test_group_column_optional(tmp_path: Path) -> None:
    p = tmp_path / "s.tsv"
    p.write_text(
        "sample_id\tbam\tstrandedness\na\ta.bam\trf\nb\tb.bam\trf\n",
        encoding="utf-8",
    )
    rows = load_samples(p)
    assert [r.group for r in rows] == ["", ""]
    assert [r.replicate for r in rows] == [1, 2]


def test_treat_only_sheet(tmp_path: Path) -> None:
    from txnova.samples import can_run_de, has_contrast

    p = tmp_path / "s.tsv"
    p.write_text(
        "sample_id\tbam\tgroup\tstrandedness\nt1\tt1.bam\ttreat\trf\nt2\tt2.bam\ttreat\trf\n",
        encoding="utf-8",
    )
    rows = load_samples(p)
    assert has_contrast(rows) is False
    assert can_run_de(rows) is False


def test_auto_replicate(tmp_path: Path) -> None:
    p = tmp_path / "s.tsv"
    p.write_text(
        "sample_id\tbam\tgroup\tstrandedness\n"
        "c1\tc1.bam\tcontrol\trf\n"
        "c2\tc2.bam\tcontrol\trf\n"
        "t1\tt1.bam\ttreat\trf\n",
        encoding="utf-8",
    )
    rows = load_samples(p)
    assert [r.replicate for r in rows] == [1, 2, 1]
