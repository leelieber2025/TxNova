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
