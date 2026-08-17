"""Evidence for loci with no gene name: conservation + Pfam.

Only `named_overlap=none`. UCSC phyloP / phastCons on exons (PhyloCSF is
not on the mm39 public track list). HMMER hmmscan against Pfam for the
longest ORF, if any. Network failures are recorded; they do not fail the run.
"""

from __future__ import annotations

import json
import re
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any, Callable, Iterable

import pandas as pd

from txnova.fold import parse_peptides_fa
from txnova.logging import get_logger

log = get_logger("txnova.orphan")

UCSC_TRACK = "https://api.genome.ucsc.edu/getData/track"
HMMER_SEARCH = "https://www.ebi.ac.uk/Tools/hmmer/api/v1/search/hmmscan"
HMMER_RESULT = "https://www.ebi.ac.uk/Tools/hmmer/api/v1/result/{tid}"

# HMMER inclusion threshold
MAX_EVALUE = 0.01
TOP_N = 5

NONE_VALUES = {"", "none", "nan", "na"}
_PFAM_VER = re.compile(r"^(PF\d+)\.\d+$", re.I)

TRACKS = {
    "mm39": {"phylop": "phyloP35way", "phastcons": "phastCons35way"},
}

ORPHAN_COLUMNS = [
    "locus_id",
    "chrom",
    "start",
    "end",
    "strand",
    "n_exons",
    "length_nt",
    "named_overlap",
    "longest_orf_aa",
    "coding_label",
    "phylop_n",
    "phylop_mean",
    "phylop_max",
    "phylop_frac_pos",
    "phastcons_n",
    "phastcons_mean",
    "pfam_name",
    "pfam_acc",
    "pfam_evalue",
    "pfam_score",
    "pfam_desc",
    "n_pfam",
    "error",
]

FetchFn = Callable[[str, str, str, int, int], list[dict]]
ScanFn = Callable[..., list[dict]]


def ucsc_genome(assembly: str) -> str | None:
    a = (assembly or "").strip()
    if a in {"GRCm39", "mm39"}:
        return "mm39"
    return None


def is_naming_none(value: Any) -> bool:
    if value is None:
        return True
    if isinstance(value, float) and pd.isna(value):
        return True
    return str(value).strip().lower() in NONE_VALUES


def parse_exons(exon_structure: Any) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    for part in str(exon_structure or "").split(","):
        part = part.strip()
        if not part or part.lower() in NONE_VALUES or part == ".":
            continue
        if "-" not in part:
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


def _num(value: Any) -> float | None:
    if value is None or value == "":
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def summarize_wig(items: Iterable[dict]) -> dict[str, Any]:
    vals: list[float] = []
    for it in items:
        try:
            start = int(it["start"])
            end = int(it["end"])
            value = float(it["value"])
        except (KeyError, TypeError, ValueError):
            continue
        n = max(0, end - start)
        if n <= 1:
            vals.append(value)
        else:
            vals.extend([value] * n)
    if not vals:
        return {"n": 0, "mean": None, "max": None, "frac_pos": None}
    return {
        "n": len(vals),
        "mean": sum(vals) / len(vals),
        "max": max(vals),
        "frac_pos": sum(1 for v in vals if v > 0) / len(vals),
    }


def _http_json(url: str, *, data: bytes | None = None, timeout: int = 60) -> dict:
    headers = {"User-Agent": "txnova/0.4", "Accept": "application/json"}
    if data is not None:
        headers["Content-Type"] = "application/json"
    req = urllib.request.Request(url, data=data, method="POST" if data else "GET", headers=headers)
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode())


def fetch_track(genome: str, track: str, chrom: str, start_0: int, end_0: int) -> list[dict]:
    url = (
        f"{UCSC_TRACK}?genome={genome};track={track};chrom={chrom}"
        f";start={start_0};end={end_0};maxItemsOutput=1000000"
    )
    payload = _http_json(url, timeout=90)
    raw = payload.get(track)
    if isinstance(raw, dict):
        raw = raw.get(chrom) or next(iter(raw.values()), [])
    if not isinstance(raw, list):
        raise RuntimeError(f"UCSC {track} empty for {chrom}:{start_0}-{end_0}")
    return raw


def score_conservation(
    chrom: str,
    exons: list[tuple[int, int]],
    *,
    genome: str | None,
    fetch: FetchFn = fetch_track,
) -> dict[str, Any]:
    empty = {
        "phylop_n": 0,
        "phylop_mean": None,
        "phylop_max": None,
        "phylop_frac_pos": None,
        "phastcons_n": 0,
        "phastcons_mean": None,
        "error": "",
    }
    tracks = TRACKS.get(genome or "")
    if not tracks:
        empty["error"] = f"no UCSC conservation track for {genome or 'unknown'}"
        return empty
    if not exons:
        empty["error"] = "no exons"
        return empty
    phy_items: list[dict] = []
    ph_items: list[dict] = []
    errors: list[str] = []
    for start, end in exons:
        start_0, end_0 = start - 1, end
        try:
            phy_items.extend(fetch(genome, tracks["phylop"], chrom, start_0, end_0))
        except (
            urllib.error.URLError,
            TimeoutError,
            OSError,
            RuntimeError,
            json.JSONDecodeError,
            KeyError,
        ) as e:
            errors.append(f"phyloP {chrom}:{start}-{end}: {e}")
        try:
            ph_items.extend(fetch(genome, tracks["phastcons"], chrom, start_0, end_0))
        except (
            urllib.error.URLError,
            TimeoutError,
            OSError,
            RuntimeError,
            json.JSONDecodeError,
            KeyError,
        ) as e:
            errors.append(f"phastCons {chrom}:{start}-{end}: {e}")
    phy = summarize_wig(phy_items)
    ph = summarize_wig(ph_items)
    return {
        "phylop_n": phy["n"],
        "phylop_mean": phy["mean"],
        "phylop_max": phy["max"],
        "phylop_frac_pos": phy["frac_pos"],
        "phastcons_n": ph["n"],
        "phastcons_mean": ph["mean"],
        "error": "; ".join(errors),
    }


def parse_hmmscan_hits(payload: dict) -> list[dict]:
    result = payload.get("result") if isinstance(payload.get("result"), dict) else payload
    hits = (result or {}).get("hits") or payload.get("hits") or []
    out: list[dict] = []
    for h in hits:
        if not isinstance(h, dict):
            continue
        ev = _num(h.get("evalue") if h.get("evalue") is not None else h.get("E"))
        if ev is None or ev > MAX_EVALUE:
            continue
        meta = h.get("metadata") if isinstance(h.get("metadata"), dict) else {}
        acc = str(meta.get("accession") or h.get("acc") or h.get("accession") or "")
        m = _PFAM_VER.match(acc)
        if m:
            acc = m.group(1)
        out.append(
            {
                "name": str(meta.get("identifier") or h.get("name") or h.get("hmmname") or ""),
                "acc": acc,
                "evalue": ev,
                "score": _num(h.get("score") if h.get("score") is not None else h.get("bits")),
                "desc": str(meta.get("description") or h.get("desc") or h.get("description") or ""),
            }
        )
    out.sort(key=lambda r: r["evalue"])
    return out


def hmmscan(seq: str, *, name: str = "query", timeout_s: int = 180) -> list[dict]:
    seq = "".join(str(seq).split())
    if not seq:
        return []
    ticket = _http_json(
        HMMER_SEARCH,
        data=json.dumps({"input": f">{name}\n{seq}", "database": "pfam"}).encode(),
        timeout=60,
    )
    tid = ticket.get("id")
    if not tid:
        raise RuntimeError(f"no HMMER ticket: {ticket}")
    deadline = time.time() + timeout_s
    while time.time() < deadline:
        payload = _http_json(HMMER_RESULT.format(tid=tid), timeout=45)
        status = str(payload.get("status") or "").upper()
        if status in {"SUCCESS", "COMPLETE", "COMPLETED", "DONE"}:
            return parse_hmmscan_hits(payload)
        if status in {"ERROR", "FAILED", "FAILURE"}:
            raise RuntimeError(f"HMMER failed: {status}")
        time.sleep(3)
    raise TimeoutError(f"HMMER timed out for {name}")


def collect_orphan_rows(frames: Iterable[pd.DataFrame]) -> pd.DataFrame:
    keep: list[pd.DataFrame] = []
    for df in frames:
        if df is None or df.empty:
            continue
        if "named_overlap" in df.columns:
            keep.append(df.loc[df["named_overlap"].map(is_naming_none)].copy())
        else:
            keep.append(df.copy())
    if not keep:
        return pd.DataFrame()
    out = pd.concat(keep, ignore_index=True)
    if "locus_id" in out.columns:
        out = out.drop_duplicates("locus_id")
    return out


def _empty_table() -> pd.DataFrame:
    return pd.DataFrame(columns=ORPHAN_COLUMNS)


def annotate_orphans(
    tables: list[Path],
    peptides_fa: Path | None,
    dest: Path,
    *,
    assembly: str = "GRCm39",
    min_orf_aa: int = 50,
    fetch: FetchFn | None = None,
    scan: ScanFn | None = None,
    only: Iterable[str] | None = None,
) -> pd.DataFrame:
    frames = [pd.read_csv(p, sep="\t") for p in tables if p.is_file()]
    src = collect_orphan_rows(frames)
    if only is not None and not src.empty and "locus_id" in src.columns:
        wanted = {str(x) for x in only}
        src = src.loc[src["locus_id"].astype(str).isin(wanted)].copy()
    if not src.empty and "longest_orf_aa" in src.columns:
        aa = pd.to_numeric(src["longest_orf_aa"], errors="coerce").fillna(0)
        src = src.loc[aa >= min_orf_aa].copy()
    dest.parent.mkdir(parents=True, exist_ok=True)
    if src.empty:
        out = _empty_table()
        out.to_csv(dest, sep="\t", index=False)
        return out
    peptides = parse_peptides_fa(peptides_fa) if peptides_fa is not None else {}
    genome = ucsc_genome(assembly)
    fetch_fn = fetch or fetch_track
    scan_fn = scan or hmmscan
    rows: list[dict] = []
    for rec in src.to_dict(orient="records"):
        locus = str(rec.get("locus_id") or "")
        chrom = str(rec.get("chrom") or "")
        cons = score_conservation(
            chrom,
            parse_exons(rec.get("exon_structure")),
            genome=genome,
            fetch=fetch_fn,
        )
        pep = peptides.get(locus, "")
        hits: list[dict] = []
        hmm_err = ""
        if pep:
            try:
                hits = list(scan_fn(pep, name=locus)[:TOP_N])
                log.info("%s hmmscan %s hits", locus, len(hits))
            except Exception as e:
                hmm_err = str(e)
                log.warning("hmmscan failed for %s: %s", locus, e)
        top = hits[0] if hits else {}
        err = "; ".join(x for x in (cons.get("error") or "", hmm_err) if x)
        rows.append(
            {
                "locus_id": locus,
                "chrom": chrom,
                "start": rec.get("start"),
                "end": rec.get("end"),
                "strand": rec.get("strand"),
                "n_exons": rec.get("n_exons"),
                "length_nt": rec.get("length_nt"),
                "named_overlap": rec.get("named_overlap")
                if not is_naming_none(rec.get("named_overlap"))
                else "none",
                "longest_orf_aa": rec.get("longest_orf_aa"),
                "coding_label": rec.get("coding_label"),
                "phylop_n": cons["phylop_n"],
                "phylop_mean": cons["phylop_mean"],
                "phylop_max": cons["phylop_max"],
                "phylop_frac_pos": cons["phylop_frac_pos"],
                "phastcons_n": cons["phastcons_n"],
                "phastcons_mean": cons["phastcons_mean"],
                "pfam_name": top.get("name") or "",
                "pfam_acc": top.get("acc") or "",
                "pfam_evalue": top.get("evalue") if top else "",
                "pfam_score": top.get("score") if top else "",
                "pfam_desc": top.get("desc") or "",
                "n_pfam": len(hits),
                "error": err,
            }
        )
        log.info(
            "%s phyloP_mean=%s pfam=%s",
            locus,
            cons["phylop_mean"],
            top.get("name") or "",
        )
    out = pd.DataFrame(rows)
    for col in ORPHAN_COLUMNS:
        if col not in out.columns:
            out[col] = ""
    out = out[ORPHAN_COLUMNS]
    out.to_csv(dest, sep="\t", index=False)
    return out
