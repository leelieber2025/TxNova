"""Join residual splice leak rows to candidate / unnamed tables."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from txnova.logging import get_logger

log = get_logger("txnova.leak")

LEAK_FILENAME = "leak.tsv"


def _locus_set(path: Path) -> set[str]:
    if not path.is_file():
        return set()
    df = pd.read_csv(path, sep="\t")
    if df.empty or "locus_id" not in df.columns:
        return set()
    return set(df["locus_id"].astype(str))


def annotate_leak(
    leak_tsv: Path,
    *,
    gates: Path,
    unnamed: Path,
    candidates: Path,
) -> pd.DataFrame:
    if not leak_tsv.is_file() or leak_tsv.stat().st_size == 0:
        return pd.DataFrame()
    df = pd.read_csv(leak_tsv, sep="\t")
    if df.empty:
        return df
    g = _locus_set(gates)
    u = _locus_set(unnamed)
    c = _locus_set(candidates)
    loc = df["merged_locus"].fillna("").astype(str) if "merged_locus" in df.columns else ""
    df["in_gates"] = loc.map(lambda x: x in g and x not in {"", "nan", "NA"})
    df["in_unnamed"] = loc.map(lambda x: x in u and x not in {"", "nan", "NA"})
    df["in_candidates"] = loc.map(lambda x: x in c and x not in {"", "nan", "NA"})
    df.to_csv(leak_tsv, sep="\t", index=False)
    n_u = int((df["status"] == "unassembled").sum()) if "status" in df.columns else 0
    n_a = int((df["status"] == "assembled_u").sum()) if "status" in df.columns else 0
    log.info("leak table: %s unassembled, %s assembled_u", n_u, n_a)
    return df
