"""Step stamps. Fingerprint = config subset + sheet sha256 + BAM size/mtime (not content)."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from txnova.logging import get_logger
from txnova.samples import SampleRow

log = get_logger("txnova.stamps")


def _sha256_text(s: str) -> str:
    return hashlib.sha256(s.encode()).hexdigest()


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def path_fingerprint(paths: list[Path]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for p in paths:
        st = p.stat() if p.is_file() else None
        out[str(p)] = {
            "size": None if st is None else st.st_size,
            "mtime": None if st is None else int(st.st_mtime),
        }
    return out


def bam_fingerprint(rows: list[SampleRow]) -> dict[str, Any]:
    out = {}
    for r in rows:
        st = r.bam.stat() if r.bam.is_file() else None
        out[r.sample_id] = {
            "path": str(r.bam),
            "size": None if st is None else st.st_size,
            "mtime": None if st is None else int(st.st_mtime),
        }
    return out


def write_stamp(
    path: Path,
    fingerprint: dict[str, Any],
    *,
    outputs: dict[str, Any] | None = None,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload: dict[str, Any] = {"fingerprint": fingerprint}
    if outputs is not None:
        payload["outputs"] = outputs
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def stamp_matches(path: Path, expected: dict[str, Any]) -> bool:
    if not path.is_file():
        return False
    try:
        got = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return False
    if not isinstance(got, dict) or "fingerprint" not in got:
        return False
    return got["fingerprint"] == expected


def stamp_outputs(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {}
    try:
        got = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {}
    if not isinstance(got, dict):
        return {}
    raw = got.get("outputs")
    return raw if isinstance(raw, dict) else {}


def clear_stamps(outdir: Path) -> None:
    d = outdir / "stamps"
    if not d.is_dir():
        return
    for p in d.glob("*.json"):
        p.unlink()
        log.info("removed stamp %s", p.name)


def sheet_sha(path: Path) -> str:
    return _sha256_file(path) if path.is_file() else ""


def config_subset_hash(obj: Any) -> str:
    return _sha256_text(json.dumps(obj, sort_keys=True, default=str))
