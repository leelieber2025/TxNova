"""Join residual splice leak rows to candidate / unnamed / shared tables."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from txnova.logging import get_logger

log = get_logger("txnova.leak")

LEAK_FILENAME = "leak.tsv"
SHARED_FILENAME = "candidates.shared.tsv"


def _locus_set(path: Path) -> set[str]:
    if not path.is_file():
        return set()
    df = pd.read_csv(path, sep="\t")
    if df.empty or "locus_id" not in df.columns:
        return set()
    return set(df["locus_id"].astype(str))


def _read_table(path: Path) -> pd.DataFrame:
    if not path.is_file() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        df = pd.read_csv(path, sep="\t")
    except pd.errors.EmptyDataError:
        return pd.DataFrame()
    return df if df is not None else pd.DataFrame()


def write_shared(
    *,
    residual: Path,
    structure_tables: list[Path],
    dest: Path,
) -> pd.DataFrame:
    """Structure-pass loci whose harvest junctions are also in control."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    res = _read_table(residual)
    frames: list[pd.DataFrame] = []
    for p in structure_tables:
        df = _read_table(p)
        if df.empty or "locus_id" not in df.columns:
            continue
        frames.append(df)
    if not frames:
        dest.write_text("locus_id\n", encoding="utf-8")
        return pd.DataFrame()
    src = pd.concat(frames, ignore_index=True)
    src = src.drop_duplicates("locus_id", keep="first")
    if res.empty or "residual_id" not in res.columns:
        dest.write_text("\t".join(src.columns.astype(str)) + "\n", encoding="utf-8")
        return pd.DataFrame(columns=list(src.columns))
    ctrl = pd.to_numeric(res.get("control_max"), errors="coerce").fillna(0)
    shared_res = res.loc[ctrl > 0].copy()
    if shared_res.empty:
        dest.write_text("\t".join(src.columns.astype(str)) + "\n", encoding="utf-8")
        return pd.DataFrame(columns=list(src.columns))
    ids = set(shared_res["residual_id"].astype(str))
    out = src.loc[src["locus_id"].astype(str).isin(ids)].copy()
    extra = shared_res.set_index(shared_res["residual_id"].astype(str))
    keep = [c for c in ("control_max", "control_n_detected", "cohort") if c in extra.columns]
    if keep:
        mapped = extra[keep]
        for col in keep:
            out[f"residual_{col}"] = out["locus_id"].astype(str).map(mapped[col])
    dest.parent.mkdir(parents=True, exist_ok=True)
    if out.empty:
        dest.write_text("\t".join(src.columns.astype(str)) + "\n", encoding="utf-8")
        log.info("shared splice table: 0 loci")
        return out
    out.to_csv(dest, sep="\t", index=False)
    log.info("shared splice table: %s structure-pass loci also in control", len(out))
    return out


def exclude_shared_from_finals(final_path: Path, shared_path: Path, rows) -> int:
    """Treatment-specific finals cannot keep a shared-harvest locus."""
    from txnova.gates import write_candidates

    if not final_path.is_file():
        return 0
    cand = _read_table(final_path)
    if cand.empty or "locus_id" not in cand.columns:
        return 0
    drop = _locus_set(shared_path)
    if not drop:
        return 0
    keep = ~cand["locus_id"].astype(str).isin(drop)
    n_drop = int((~keep).sum())
    if n_drop:
        write_candidates(cand.loc[keep].copy(), final_path, rows)
        log.info("finals: dropped %s shared-harvest loci", n_drop)
    return n_drop


def annotate_leak(
    leak_tsv: Path,
    *,
    gates: Path,
    unnamed: Path,
    candidates: Path,
    shared: Path | None = None,
) -> pd.DataFrame:
    if not leak_tsv.is_file() or leak_tsv.stat().st_size == 0:
        return pd.DataFrame()
    df = pd.read_csv(leak_tsv, sep="\t")
    if df.empty:
        return df
    g = _locus_set(gates)
    u = _locus_set(unnamed)
    c = _locus_set(candidates)
    s = _locus_set(shared) if shared is not None else set()
    loc = df["merged_locus"].fillna("").astype(str) if "merged_locus" in df.columns else ""
    df["in_gates"] = loc.map(lambda x: x in g and x not in {"", "nan", "NA"})
    df["in_unnamed"] = loc.map(lambda x: x in u and x not in {"", "nan", "NA"})
    df["in_candidates"] = loc.map(lambda x: x in c and x not in {"", "nan", "NA"})
    df["in_shared"] = loc.map(lambda x: x in s and x not in {"", "nan", "NA"})
    df.to_csv(leak_tsv, sep="\t", index=False)
    n_u = int((df["status"] == "unassembled").sum()) if "status" in df.columns else 0
    n_a = int((df["status"] == "assembled_u").sum()) if "status" in df.columns else 0
    n_sh = int((df["cohort"] == "shared").sum()) if "cohort" in df.columns else 0
    log.info("leak table: %s unassembled, %s assembled_u, %s shared", n_u, n_a, n_sh)
    return df
