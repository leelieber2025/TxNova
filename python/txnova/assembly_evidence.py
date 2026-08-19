"""Layer 1: which sample StringTie GTFs proposed a merge locus."""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from txnova.gates import write_candidates
from txnova.logging import get_logger
from txnova.samples import SampleRow

log = get_logger("txnova.assembly_evidence")

_ATTR = re.compile(r'(\S+)\s+"([^"]*)"')


def parse_stringtie_transcripts(path: Path) -> list[dict]:
    rows: list[dict] = []
    if not path.is_file():
        return rows
    with path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "transcript":
                continue
            attrs = dict(_ATTR.findall(cols[8]))
            try:
                cov = float(attrs["cov"]) if attrs.get("cov") else None
            except ValueError:
                cov = None
            rows.append(
                {
                    "chrom": cols[0],
                    "start": int(cols[3]),
                    "end": int(cols[4]),
                    "strand": cols[6],
                    "cov": cov,
                }
            )
    return rows


def _overlap(a0: int, a1: int, b0: int, b1: int) -> bool:
    return not (a1 < b0 or b1 < a0)


def _hits_for_locus(
    txs: list[dict],
    chrom: str,
    start: int,
    end: int,
    strand: str,
    unstranded: bool,
) -> list[dict]:
    out = []
    for t in txs:
        if t["chrom"] != chrom:
            continue
        if not unstranded and strand not in {".", ""} and t["strand"] not in {".", strand}:
            continue
        if _overlap(start, end, t["start"], t["end"]):
            out.append(t)
    return out


def attach_assembly_evidence(
    cand_path: Path,
    rows: list[SampleRow],
    assembly_dir: Path,
    *,
    unstranded: bool = False,
) -> None:
    if not cand_path.is_file():
        return
    cand = pd.read_csv(cand_path, sep="\t")
    extra_cols = []
    for r in rows:
        extra_cols.extend(
            [
                f"{r.sample_id}_assembled",
                f"{r.sample_id}_asm_n_isoforms",
                f"{r.sample_id}_asm_max_cov",
                f"{r.sample_id}_asm_span",
            ]
        )
    extra_cols.extend(
        ["assembled_in", "n_assembled_control", "n_assembled_treat", "asm_span_max_delta_bp"]
    )
    if cand.empty:
        for c in extra_cols:
            if c not in cand.columns:
                cand[c] = pd.Series(dtype=object)
        write_candidates(cand, cand_path, rows)
        return
    by_sample = {
        r.sample_id: parse_stringtie_transcripts(assembly_dir / f"{r.sample_id}.gtf") for r in rows
    }
    controls = {r.sample_id for r in rows if r.group == "control"}
    treats = {r.sample_id for r in rows if r.group == "treat"}

    extra: list[dict] = []
    for rec in cand.to_dict(orient="records"):
        chrom = str(rec["chrom"])
        start = int(rec["start"])
        end = int(rec["end"])
        strand = str(rec["strand"])
        assembled: list[str] = []
        n_ctrl = 0
        n_treat = 0
        max_delta = 0
        row: dict = {}
        for r in rows:
            hits = _hits_for_locus(by_sample[r.sample_id], chrom, start, end, strand, unstranded)
            proposed = bool(hits)
            row[f"{r.sample_id}_assembled"] = proposed
            row[f"{r.sample_id}_asm_n_isoforms"] = len(hits)
            covs = [h["cov"] for h in hits if h["cov"] is not None]
            row[f"{r.sample_id}_asm_max_cov"] = max(covs) if covs else ""
            if hits:
                s0 = min(h["start"] for h in hits)
                s1 = max(h["end"] for h in hits)
                row[f"{r.sample_id}_asm_span"] = f"{s0}-{s1}"
                delta = max(abs(s0 - start), abs(s1 - end))
                if delta > max_delta:
                    max_delta = delta
                assembled.append(r.sample_id)
                if r.sample_id in controls:
                    n_ctrl += 1
                if r.sample_id in treats:
                    n_treat += 1
            else:
                row[f"{r.sample_id}_asm_span"] = ""
        row["assembled_in"] = ",".join(assembled)
        row["n_assembled_control"] = n_ctrl
        row["n_assembled_treat"] = n_treat
        row["asm_span_max_delta_bp"] = max_delta if assembled else ""
        extra.append(row)

    out = pd.concat([cand.reset_index(drop=True), pd.DataFrame(extra)], axis=1)
    write_candidates(out, cand_path, rows)
    log.info("assembly evidence on %s candidates", len(out))
