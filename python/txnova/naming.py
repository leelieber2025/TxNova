"""Label loci with a comprehensive gene GTF (naming, not class u)."""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path

import pandas as pd

from txnova.errors import TxNovaError
from txnova.gates import NAMING_COLUMNS, write_candidates
from txnova.samples import SampleRow

_ATTR = re.compile(r"(\S+)\s+\"([^\"]*)\"")

NAMING_NONE = "none"
NAMING_PROTEIN = "protein_coding"
NAMING_ANNOTATED = "annotated"


def chrom_key(name: str) -> str:
    n = str(name)
    n = n.removeprefix("chr")
    if n.upper() in {"M", "MT"}:
        return "MT"
    return n


def load_gene_bodies(gtf: Path) -> dict[str, list[tuple[int, int, str, str, str, str]]]:
    if not gtf.is_file():
        raise TxNovaError(f"naming annotation not found: {gtf}")
    by_chrom: dict[str, list[tuple[int, int, str, str, str, str]]] = defaultdict(list)
    # gene_id → (chrom_key, start, end, strand, name, type, id)
    tx_span: dict[str, tuple[str, int, int, str, str, str, str]] = {}
    with gtf.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] not in {"gene", "transcript", "exon"}:
                continue
            try:
                start, end = int(p[3]), int(p[4])
            except ValueError:
                continue
            attrs = dict(_ATTR.findall(p[8]))
            rec = (
                start,
                end,
                p[6],
                attrs.get("gene_name", ""),
                attrs.get("gene_type") or attrs.get("gene_biotype", ""),
                attrs.get("gene_id", ""),
            )
            ck = chrom_key(p[0])
            if p[2] == "gene":
                by_chrom[ck].append(rec)
                continue
            gid = rec[5]
            if not gid:
                continue
            prev = tx_span.get(gid)
            if prev is None:
                tx_span[gid] = (ck, start, end, rec[2], rec[3], rec[4], gid)
            else:
                tx_span[gid] = (
                    prev[0],
                    min(prev[1], start),
                    max(prev[2], end),
                    prev[3],
                    prev[4] or rec[3],
                    prev[5] or rec[4],
                    gid,
                )
    if by_chrom:
        return by_chrom
    for ck, start, end, strand, name, typ, gid in tx_span.values():
        by_chrom[ck].append((start, end, strand, name, typ, gid))
    return by_chrom


def _overlap(a0: int, a1: int, b0: int, b1: int) -> bool:
    return a0 <= b1 and b0 <= a1


def classify_locus(
    chrom: str,
    start: int,
    end: int,
    strand: str,
    genes: dict[str, list[tuple[int, int, str, str, str, str]]],
) -> dict[str, str]:
    empty = {
        "named_gene_name": "",
        "named_gene_id": "",
        "named_gene_type": "",
        "named_overlap": NAMING_NONE,
    }
    hits = [g for g in genes.get(chrom_key(chrom), []) if _overlap(start, end, g[0], g[1])]
    if not hits:
        return empty

    def pick(pred) -> tuple[int, int, str, str, str, str] | None:
        matched = [g for g in hits if pred(g)]
        if not matched:
            return None
        return min(matched, key=lambda g: (g[0], g[5]))

    same_pc = pick(lambda g: g[2] == strand and g[4] == "protein_coding")
    if same_pc:
        return {
            "named_gene_name": same_pc[3],
            "named_gene_id": same_pc[5],
            "named_gene_type": same_pc[4],
            "named_overlap": NAMING_PROTEIN,
        }
    any_hit = pick(lambda g: g[2] == strand) or pick(lambda _g: True)
    if any_hit is None:
        return empty
    return {
        "named_gene_name": any_hit[3],
        "named_gene_id": any_hit[5],
        "named_gene_type": any_hit[4],
        "named_overlap": NAMING_ANNOTATED,
    }


def annotate_table(df: pd.DataFrame, genes: dict) -> pd.DataFrame:
    out = df.copy()
    for c in NAMING_COLUMNS:
        if c in out.columns:
            out = out.drop(columns=[c])
    if out.empty or "chrom" not in out.columns:
        for c in NAMING_COLUMNS:
            out[c] = "" if c != "named_overlap" else NAMING_NONE
        return out
    rows = []
    for rec in out.itertuples(index=False):
        try:
            start, end = int(rec.start), int(rec.end)
        except (TypeError, ValueError):
            rows.append(
                {
                    "named_gene_name": "",
                    "named_gene_id": "",
                    "named_gene_type": "",
                    "named_overlap": NAMING_NONE,
                }
            )
            continue
        rows.append(classify_locus(str(rec.chrom), start, end, str(rec.strand), genes))
    labels = pd.DataFrame(rows)
    return pd.concat([out.reset_index(drop=True), labels], axis=1)


def label_candidate_file(
    path: Path,
    rows: list[SampleRow],
    genes: dict[str, list[tuple[int, int, str, str, str, str]]],
) -> None:
    df = pd.read_csv(path, sep="\t") if path.is_file() else pd.DataFrame()
    write_candidates(annotate_table(df, genes), path, rows)
