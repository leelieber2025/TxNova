from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
from pydantic import BaseModel, ConfigDict, Field

from txnova.errors import TxNovaError

SAMPLE_ID_RE = re.compile(r"^[A-Za-z0-9._-]+$")
REQUIRED_COLS = ("sample_id", "bam", "group", "strandedness")


class SampleRow(BaseModel):
    model_config = ConfigDict(extra="forbid")
    sample_id: str
    bam: Path
    group: str
    strandedness: str
    replicate: int = Field(ge=1)


def load_samples(path: Path) -> list[SampleRow]:
    path = path.resolve()
    if not path.is_file():
        raise TxNovaError(f"sample sheet not found: {path}")
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise TxNovaError(f"sample sheet missing columns: {missing}")
    extra = [c for c in df.columns if c not in (*REQUIRED_COLS, "replicate")]
    if extra:
        raise TxNovaError(f"sample sheet has unknown columns {extra}; extra=forbid")
    if df.empty:
        raise TxNovaError("sample sheet has no rows")

    rows: list[SampleRow] = []
    seen: set[str] = set()
    group_ord: dict[str, int] = {"control": 0, "treat": 0}
    for i, rec in enumerate(df.to_dict(orient="records"), start=2):
        sid = str(rec["sample_id"]).strip()
        if not SAMPLE_ID_RE.match(sid):
            raise TxNovaError(f"line {i}: invalid sample_id {sid!r}")
        if sid in seen:
            raise TxNovaError(f"duplicate sample_id {sid}")
        seen.add(sid)
        group = str(rec["group"]).strip()
        if group not in {"control", "treat"}:
            raise TxNovaError(
                f"line {i}: group must be control|treat, got {group!r} (no case/ko/wt aliases)"
            )
        strand = str(rec["strandedness"]).strip()
        if strand not in {"unstranded", "fr", "rf"}:
            raise TxNovaError(f"line {i}: strandedness must be unstranded|fr|rf, got {strand!r}")
        bam_raw = str(rec["bam"]).strip()
        bam = Path(bam_raw)
        if not bam.is_absolute():
            bam = (path.parent / bam).resolve()
        rep_raw = str(rec.get("replicate", "")).strip()
        if rep_raw:
            try:
                replicate = int(rep_raw)
            except ValueError as e:
                raise TxNovaError(f"line {i}: replicate must be an integer") from e
        else:
            group_ord[group] += 1
            replicate = group_ord[group]
        rows.append(
            SampleRow(
                sample_id=sid,
                bam=bam,
                group=group,
                strandedness=strand,
                replicate=replicate,
            )
        )
    return rows


def samples_to_jsonable(rows: list[SampleRow]) -> list[dict]:
    return [
        {
            "sample_id": r.sample_id,
            "bam": str(r.bam),
            "group": r.group,
            "strandedness": r.strandedness,
            "replicate": r.replicate,
        }
        for r in rows
    ]
