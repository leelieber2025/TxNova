"""RepeatMasker overlap on residual exons. Used as an optional detection gate."""

from __future__ import annotations

from pathlib import Path

from txnova.errors import TxNovaError


def _chrom_key(name: str) -> str:
    n = str(name)
    n = n.removeprefix("chr")
    if n.upper() in {"M", "MT"}:
        return "MT"
    return n


def parse_exons(exon_structure: str) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    for part in str(exon_structure or "").split(","):
        part = part.strip()
        if not part or "-" not in part:
            continue
        a, b = part.split("-", 1)
        try:
            start, end = int(a), int(b)
        except ValueError:
            continue
        if end < start:
            start, end = end, start
        out.append((start, end))
    return out


def load_rmsk(path: Path) -> dict[str, list[tuple[int, int, str]]]:
    path = Path(path)
    if not path.is_file():
        raise TxNovaError(f"RepeatMasker BED not found: {path}")
    by: dict[str, list[tuple[int, int, str]]] = {}
    try:
        fh = path.open()
    except OSError as e:
        raise TxNovaError(f"RepeatMasker BED not found: {path}") from e
    with fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            chrom, start, end = p[0], int(p[1]), int(p[2])
            fam = p[4] if len(p) > 4 else (p[3] if len(p) > 3 else "")
            rec = (start, end, fam)
            ck = _chrom_key(chrom)
            by.setdefault(ck, []).append(rec)
            if chrom != ck:
                by.setdefault(chrom, []).append(rec)
    for chrom in by:
        by[chrom].sort()
    return by


def _overlap(a0: int, a1: int, b0: int, b1: int) -> int:
    return max(0, min(a1, b1) - max(a0, b0))


def rmsk_frac(
    chrom: str,
    exon_structure: str,
    length_nt: int,
    idx: dict[str, list[tuple[int, int, str]]],
) -> float:
    if length_nt <= 0:
        return 0.0
    intervals = idx.get(_chrom_key(chrom), []) or idx.get(chrom, [])
    if not intervals:
        return 0.0
    total = 0
    for es, ee in parse_exons(exon_structure):
        e0, e1 = es - 1, ee
        for rs, re, _fam in intervals:
            if re <= e0:
                continue
            if rs >= e1:
                break
            total += _overlap(e0, e1, rs, re)
    return total / length_nt
