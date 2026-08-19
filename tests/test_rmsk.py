from pathlib import Path

from txnova.rmsk import load_rmsk, rmsk_frac


def test_rmsk_frac(tmp_path: Path) -> None:
    bed = tmp_path / "rmsk.bed"
    bed.write_text("chr1\t1000\t1100\tAlu\tAlu\n")
    idx = load_rmsk(bed)
    # one 200 nt exon, 100 nt overlap
    assert abs(rmsk_frac("chr1", "1001-1200", 200, idx) - 0.5) < 1e-9
    assert rmsk_frac("chr2", "1001-1200", 200, idx) == 0.0
