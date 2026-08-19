from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from txnova.config import TxNovaConfig
from txnova.errors import TxNovaError
from txnova.logging import get_logger
from txnova.samples import SampleRow, can_run_de, samples_to_jsonable
from txnova.staging import StagingDir

log = get_logger("txnova.preflight")


def _core():
    try:
        from txnova import _core as core
    except ImportError as e:
        raise TxNovaError(
            "txnova._core is not built. Install with `pip install -e .` (maturin compiles the Rust engine)."
        ) from e
    return core


def build_samples_json(cfg: TxNovaConfig, rows: list[SampleRow]) -> str:
    payload = {
        "assembly": cfg.genome.assembly,
        "library_layout": cfg.quantify.library_layout,
        "de_enabled": cfg.de.enabled,
        "treat_min_detected_replicates": cfg.filters.treat_min_detected_replicates,
        "threads": cfg.threads,
        "samples": samples_to_jsonable(rows),
    }
    return json.dumps(payload)


def run_preflight(cfg: TxNovaConfig, rows: list[SampleRow]) -> dict[str, Any]:
    if cfg.de.enabled and not can_run_de(rows):
        log.info("DE skipped: need ≥2 control and ≥2 treat")
        cfg.de.enabled = False
    core = _core()
    payload = build_samples_json(cfg, rows)
    report = core.preflight_bams(payload, str(cfg.genome.fasta), str(cfg.genome.annotation))
    if not isinstance(report, dict):
        raise TxNovaError("preflight_bams returned a non-dict")
    if cfg.de.enabled:
        try:
            import pydeseq2  # noqa: F401
        except ImportError as e:
            raise TxNovaError(
                "de.enabled but pydeseq2 is not installed. pip install 'txnova'"
            ) from e
    return report


def write_preflight_json(cfg: TxNovaConfig, report: dict[str, Any]) -> Path:
    out = cfg.output_dir
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    with StagingDir(out) as st:
        st.write_text("preflight.json", text)
        st.publish(["preflight.json"])
    return out / "preflight.json"


def require_ok(report: dict[str, Any]) -> None:
    if report.get("ok"):
        return
    errors = report.get("errors") or ["preflight failed"]
    raise TxNovaError("preflight failed:\n  " + "\n  ".join(str(e) for e in errors))
