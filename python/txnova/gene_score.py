"""Rank structure-pass loci by how gene-like they are.

Cheap features follow CPAT/CPC2 (ORF length, hexamer, Fickett) plus
GENCODE-style filters (spliced, intergenic, not a TE). Penalties come
from the 6v6 screens: unplaced contigs, zero-support 'introns',
gag/pol/MLV, and intronless copies of known proteins (processed
pseudogene / NUMT).

Stage 1 scores every row without the network. Stage 2 runs phyloP +
Pfam on the top pool, applies TE/processed penalties, and re-sorts.
Only the top REPORT_N are function-scored.
"""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Any

import pandas as pd

from txnova.fold import parse_peptides_fa
from txnova.logging import get_logger
from txnova.naming import chrom_key
from txnova.orphan import hmmscan, parse_exons, score_conservation, ucsc_genome

log = get_logger("txnova.gene_score")

GENE_RANK_FILENAME = "gene_rank.tsv"
REPORT_N = 30
FUNCTION_POOL = 40


def top_rank_ids(path: Path, n: int = REPORT_N) -> list[str]:
    if not path.is_file() or path.stat().st_size == 0:
        return []
    df = pd.read_csv(path, sep="\t")
    if df.empty or "locus_id" not in df.columns:
        return []
    if "gene_rank" in df.columns:
        df = df.sort_values("gene_rank", ascending=True)
    return [str(x) for x in df["locus_id"].head(n).tolist()]


_UNPLACED = re.compile(
    r"^(GL|JH|KB|KQ|KI|MU|chrUn|Un_|NW_|NT_)",
    re.IGNORECASE,
)
_PRIMARY = re.compile(r"^chr([1-9]|1[0-9]|2[0-2]|[XY])$", re.IGNORECASE)
_TE = re.compile(
    r"gag|pol\b|rve|rnase[_ ]?h|mlvin|integrase|rvt|ltr|erv\b|transpos|dde_tnp|"
    r"retrovir|tlv_|capsid|\bcoat\b",
    re.IGNORECASE,
)
_PROCESSED = re.compile(
    r"atp.?synt|mt_atp|cox[1-3]?\b|cytochrom|nadh_dehydro|ef-hand|aif-1|"
    r"globin|histone|ubiquitin|gapdh|actin",
    re.IGNORECASE,
)


def _num(v: Any) -> float | None:
    if v is None or v == "" or (isinstance(v, float) and pd.isna(v)):
        return None
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _int(v: Any, default: int = 0) -> int:
    n = _num(v)
    if n is None:
        return default
    return int(n)


def junction_min(value: Any) -> int | None:
    """Minimum BAM splice support across listed introns. None if no list."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    s = str(value).strip()
    if s.lower() in {"", "nan", "na", "."}:
        return None
    vals: list[int] = []
    for part in s.replace(";", ",").split(","):
        part = part.strip()
        if part.isdigit() or (part.startswith("-") and part[1:].isdigit()):
            vals.append(int(part))
    if not vals:
        return None
    return min(vals)


def chrom_class(chrom: Any) -> str:
    raw = str(chrom or "")
    key = chrom_key(raw)
    probe = raw if _PRIMARY.match(raw) else f"chr{key}"
    if key.upper() != "MT" and _PRIMARY.match(probe):
        return "primary"
    if raw.lower() in {"chry"} and not _PRIMARY.match(probe):
        return "chrY"
    if _UNPLACED.match(raw) or _UNPLACED.match(key) or (raw and not raw.lower().startswith("chr")):
        return "unplaced"
    return "other"


def _te_hit(name: Any, desc: Any) -> bool:
    blob = f"{name or ''} {desc or ''}"
    return bool(_TE.search(blob))


def _processed_hit(name: Any, desc: Any) -> bool:
    blob = f"{name or ''} {desc or ''}"
    return bool(_PROCESSED.search(blob))


def stage1_score(rec: dict[str, Any]) -> tuple[float, list[str]]:
    flags: list[str] = []
    score = 0.0
    kind = chrom_class(rec.get("chrom"))
    if kind == "primary":
        score += 20
    elif kind == "chrY":
        score += 5
        flags.append("chrY")
    else:
        score -= 40
        flags.append("unplaced")

    n_ex = _int(rec.get("n_exons"), 1)
    if n_ex >= 3:
        score += 25
    elif n_ex == 2:
        score += 12
    else:
        flags.append("monoexon")

    jmin = junction_min(rec.get("junction_support"))
    if n_ex >= 2:
        if jmin is None:
            score -= 8
            flags.append("splice_unknown")
        elif jmin >= 2:
            score += 15
        elif jmin == 1:
            score += 5
            flags.append("weak_splice")
        else:
            score -= 25
            flags.append("no_splice")

    length = _int(rec.get("length_nt"), 0)
    if 200 <= length <= 4000:
        score += 8
    elif length > 8000:
        score -= 8
        flags.append("long_rna")

    start, end = _int(rec.get("start")), _int(rec.get("end"))
    span = end - start + 1 if end >= start else length
    if span > 100_000 and n_ex < 5:
        score -= 20
        flags.append("mega_span")
    elif span > 50_000 and n_ex <= 2:
        score -= 10
        flags.append("wide_span")

    aa = _int(rec.get("longest_orf_aa"), 0)
    if aa >= 50:
        score += min(20.0, aa / 20.0)
        if str(rec.get("orf_complete")).strip().lower() in {"true", "1"}:
            score += 4
    else:
        flags.append("no_orf")

    hexamer = _num(rec.get("coding_score"))
    label = str(rec.get("coding_label") or "").strip().lower()
    if label in {"coding", "hexamer_positive"} and hexamer is not None and hexamer > 0:
        score += 15
    elif hexamer is not None and hexamer > 0:
        score += 8
    elif hexamer is not None and hexamer < -0.3:
        score -= 5

    dist = _num(rec.get("nearest_distance_bp"))
    if dist is not None and dist >= 5000:
        score += min(10.0, math.log10(dist + 1))
    any_d = _num(rec.get("nearest_any_distance_bp"))
    if any_d is not None and 0 < any_d < 200:
        score -= 8
        flags.append("close_any")

    treat_n = _int(rec.get("treat_n_detected"), 0)
    if treat_n >= 5:
        score += 5
    return score, flags


def stage2_adjust(
    rec: dict[str, Any],
    *,
    pfam_name: str,
    pfam_desc: str,
    phylop_mean: float | None,
) -> tuple[float, list[str]]:
    extra = 0.0
    flags: list[str] = []
    n_ex = _int(rec.get("n_exons"), 1)
    if _te_hit(pfam_name, pfam_desc):
        extra -= 50
        flags.append("te")
    elif n_ex == 1 and _processed_hit(pfam_name, pfam_desc):
        extra -= 30
        flags.append("processed")
        if phylop_mean is not None and phylop_mean >= 0.8:
            extra -= 8
    if n_ex == 1 and phylop_mean is not None and phylop_mean >= 1.5 and pfam_name:
        if "processed" not in flags and "te" not in flags:
            extra -= 12
            flags.append("processed_like")
    return extra, flags


def _collect(tables: list[Path]) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for i, p in enumerate(tables):
        if not p.is_file() or p.stat().st_size == 0:
            continue
        df = pd.read_csv(p, sep="\t")
        if df.empty:
            continue
        df = df.copy()
        df["_src"] = i
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    out = pd.concat(frames, ignore_index=True)
    if "locus_id" in out.columns:
        out = out.sort_values("_src").drop_duplicates("locus_id", keep="first")
    return out.drop(columns=["_src"])


def rank_loci(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    scores: list[float] = []
    flags: list[str] = []
    for rec in df.to_dict(orient="records"):
        s, f = stage1_score(rec)
        scores.append(s)
        flags.append(",".join(f))
    out = df.copy()
    out["gene_score"] = scores
    out["gene_flags"] = flags
    out = out.sort_values(["gene_score", "n_exons", "length_nt"], ascending=False)
    out["gene_rank"] = range(1, len(out) + 1)
    return out


def _annotate_function(
    recs: list[dict[str, Any]],
    peptides: dict[str, str],
    *,
    assembly: str,
    scan=None,
    fetch=None,
) -> dict[str, dict[str, Any]]:
    genome = ucsc_genome(assembly)
    scan_fn = scan or hmmscan
    out: dict[str, dict[str, Any]] = {}
    for rec in recs:
        locus = str(rec.get("locus_id") or "")
        cons_kw = {"genome": genome}
        if fetch is not None:
            cons_kw["fetch"] = fetch
        cons = score_conservation(
            str(rec.get("chrom") or ""),
            parse_exons(rec.get("exon_structure")),
            **cons_kw,
        )
        pep = peptides.get(locus, "")
        hits: list[dict] = []
        err = ""
        if pep:
            try:
                hits = list(scan_fn(pep, name=locus)[:5])
            except Exception as e:
                err = str(e)
                log.warning("hmmscan failed for %s: %s", locus, e)
        top = hits[0] if hits else {}
        out[locus] = {
            "phylop_mean": cons.get("phylop_mean"),
            "phastcons_mean": cons.get("phastcons_mean"),
            "pfam_name": top.get("name") or "",
            "pfam_acc": top.get("acc") or "",
            "pfam_evalue": top.get("evalue") if top else "",
            "pfam_desc": top.get("desc") or "",
            "n_pfam": len(hits),
            "function_error": "; ".join(x for x in (str(cons.get("error") or ""), err) if x),
        }
    return out


def write_gene_rank(
    tables: list[Path],
    dest: Path,
    *,
    peptides_fa: Path | None = None,
    assembly: str = "GRCm39",
    top_n: int = REPORT_N,
    pool_n: int = FUNCTION_POOL,
    scan=None,
    fetch=None,
) -> pd.DataFrame:
    src = _collect(tables)
    dest.parent.mkdir(parents=True, exist_ok=True)
    if src.empty:
        dest.write_text(
            "gene_rank\tlocus_id\tlocus_coord\tgene_score\tgene_flags\n",
            encoding="utf-8",
        )
        return pd.DataFrame()
    ranked = rank_loci(src)
    peptides = parse_peptides_fa(peptides_fa) if peptides_fa is not None else {}
    head = ranked.head(pool_n)
    if not peptides and scan is None and fetch is None:
        cols = [
            c
            for c in (
                "gene_rank",
                "locus_id",
                "locus_coord",
                "chrom",
                "n_exons",
                "length_nt",
                "nearest_gene_name",
                "nearest_distance_bp",
                "treat_median_tpm",
                "control_max_tpm",
                "treat_n_detected",
                "junction_support",
                "longest_orf_aa",
                "coding_score",
                "coding_label",
                "gene_score",
                "gene_flags",
            )
            if c in ranked.columns
        ]
        out = ranked[cols]
        out.to_csv(dest, sep="\t", index=False)
        log.info("gene_rank: %s loci (no function pass)", len(out))
        return out
    fn = _annotate_function(
        head.to_dict(orient="records"),
        peptides,
        assembly=assembly,
        scan=scan,
        fetch=fetch,
    )
    for col in (
        "phylop_mean",
        "phastcons_mean",
        "pfam_name",
        "pfam_acc",
        "pfam_evalue",
        "pfam_desc",
        "n_pfam",
        "function_error",
    ):
        if col not in ranked.columns:
            ranked[col] = pd.Series([pd.NA] * len(ranked), dtype=object, index=ranked.index)
        else:
            ranked[col] = ranked[col].astype(object)
    extra_scores: dict[str, float] = {}
    extra_flags: dict[str, str] = {}
    for rec in head.to_dict(orient="records"):
        locus = str(rec["locus_id"])
        info = fn.get(locus) or {}
        adj, fl = stage2_adjust(
            rec,
            pfam_name=str(info.get("pfam_name") or ""),
            pfam_desc=str(info.get("pfam_desc") or ""),
            phylop_mean=_num(info.get("phylop_mean")),
        )
        extra_scores[locus] = adj
        extra_flags[locus] = ",".join(fl)
        for k, v in info.items():
            ranked.loc[ranked["locus_id"].astype(str) == locus, k] = v if v is not None else ""
    if extra_scores:
        ranked["gene_score"] = [
            float(s) + extra_scores.get(str(i), 0.0)
            for s, i in zip(ranked["gene_score"], ranked["locus_id"])
        ]
        ranked["gene_flags"] = [
            ",".join(x for x in (str(f), extra_flags.get(str(i), "")) if x)
            for f, i in zip(ranked["gene_flags"], ranked["locus_id"])
        ]
        ranked = ranked.sort_values(["gene_score", "n_exons", "length_nt"], ascending=False)
        ranked["gene_rank"] = range(1, len(ranked) + 1)
    cols = [
        c
        for c in (
            "gene_rank",
            "locus_id",
            "locus_coord",
            "chrom",
            "n_exons",
            "length_nt",
            "nearest_gene_name",
            "nearest_distance_bp",
            "treat_median_tpm",
            "control_max_tpm",
            "treat_n_detected",
            "junction_support",
            "longest_orf_aa",
            "coding_score",
            "coding_label",
            "gene_score",
            "gene_flags",
            "phylop_mean",
            "phastcons_mean",
            "pfam_name",
            "pfam_evalue",
            "pfam_desc",
            "n_pfam",
        )
        if c in ranked.columns
    ]
    out = ranked[cols]
    out.to_csv(dest, sep="\t", index=False)
    log.info(
        "gene_rank: %s loci, function on top %s, report %s", len(out), min(pool_n, len(out)), top_n
    )
    return out
