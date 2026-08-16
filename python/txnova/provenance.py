from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from txnova import __version__
from txnova.config import TxNovaConfig
from txnova.errors import TxNovaError
from txnova.samples import SampleRow
from txnova.staging import StagingDir


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _bam_meta(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise TxNovaError(f"BAM not found: {path}")
    st = path.stat()
    return {
        "path": str(path),
        "size": st.st_size,
        "mtime": int(st.st_mtime),
    }


def _core_version() -> str:
    try:
        from txnova import _core as core

        return str(core.core_version())
    except Exception:
        return "unbuilt"


def build_run_json(
    cfg: TxNovaConfig,
    rows: list[SampleRow],
    preflight: dict[str, Any],
    *,
    config_path: Path,
    phase: str = "0",
) -> dict[str, Any]:
    canonical = cfg.model_dump_canonical()
    canonical_text = json.dumps(canonical, sort_keys=True, separators=(",", ":"))
    return {
        "txnova_version": __version__,
        "txnova_core_version": _core_version(),
        "phase": phase,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "config_path": str(config_path),
        "config": canonical,
        "config_sha256": hashlib.sha256(canonical_text.encode()).hexdigest(),
        "sample_sheet": str(cfg.samples),
        "sample_sheet_sha256": _sha256_file(cfg.samples) if cfg.samples.is_file() else None,
        "genome": {
            "fasta": str(cfg.genome.fasta),
            "annotation": str(cfg.genome.annotation),
            "annotation_source": cfg.genome.annotation_source,
            "annotation_version": cfg.genome.annotation_version,
            "assembly": cfg.genome.assembly,
        },
        "bams": {r.sample_id: _bam_meta(r.bam) for r in rows},
        "sq_sha256": preflight.get("sq_sha256"),
        "aligner_family": preflight.get("aligner_family"),
        "library_layout": preflight.get("library_layout"),
        "threads": preflight.get("threads"),
        "warnings": preflight.get("warnings") or [],
        "assembler": None,
        "discovery": "not run (Phase 0)",
    }


def write_run_json(cfg: TxNovaConfig, payload: dict[str, Any]) -> Path:
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    with StagingDir(cfg.output_dir) as st:
        st.write_text("run.json", text)
        st.publish(["run.json"])
    return cfg.output_dir / "run.json"
