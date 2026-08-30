from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from txnova.config import TxNovaConfig
from txnova.errors import TxNovaError
from txnova.logging import get_logger
from txnova.naming import chrom_key
from txnova.rmsk import load_rmsk
from txnova.samples import SampleRow, can_run_de, has_contrast, n_treat, samples_to_jsonable
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


def _fai_chroms(fasta: Path) -> set[str]:
    fai = Path(str(fasta) + ".fai")
    if not fai.is_file():
        return set()
    names: set[str] = set()
    for line in fai.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.strip():
            continue
        names.add(line.split("\t")[0])
    return names


def _python_preflight_checks(
    cfg: TxNovaConfig, rows: list[SampleRow]
) -> tuple[list[str], list[str]]:
    errors: list[str] = []
    warnings: list[str] = []
    labeled = [r for r in rows if r.group in {"control", "treat"}]
    unlabeled = [r for r in rows if r.group == ""]
    if labeled and unlabeled:
        errors.append(
            "sample sheet mixes blank group with labeled control/treat; "
            "fill every group or leave every group blank"
        )
    n_t = n_treat(rows)
    if has_contrast(rows) and cfg.filters.treat_min_detected_replicates > n_t:
        errors.append(
            f"filters.treat_min_detected_replicates="
            f"{cfg.filters.treat_min_detected_replicates} exceeds n_treat={n_t}"
        )
    if cfg.genome.rmsk_bed is not None:
        rmsk = Path(cfg.genome.rmsk_bed)
        if not rmsk.is_file():
            errors.append(f"genome.rmsk_bed not found: {rmsk}")
        else:
            try:
                idx = load_rmsk(rmsk)
            except TxNovaError as e:
                errors.append(str(e))
                idx = {}
            fai_keys = {chrom_key(n) for n in _fai_chroms(cfg.genome.fasta)}
            rmsk_keys = {chrom_key(n) for n in idx}
            if fai_keys and rmsk_keys and fai_keys.isdisjoint(rmsk_keys):
                warnings.append(
                    "genome.rmsk_bed contig names are disjoint from FASTA .fai "
                    f"(BED e.g. {next(iter(idx))!r}, FASTA e.g. {next(iter(_fai_chroms(cfg.genome.fasta)))!r}); "
                    "max_rmsk_frac will be 0 for every locus"
                )
    return errors, warnings


def run_preflight(cfg: TxNovaConfig, rows: list[SampleRow]) -> dict[str, Any]:
    if cfg.de.enabled and not can_run_de(rows):
        log.info("DE skipped: need ≥2 control and ≥2 treat")
        cfg.de.enabled = False
    py_errors, py_warnings = _python_preflight_checks(cfg, rows)
    core = _core()
    payload = build_samples_json(cfg, rows)
    report = core.preflight_bams(payload, str(cfg.genome.fasta), str(cfg.genome.annotation))
    if not isinstance(report, dict):
        raise TxNovaError("preflight_bams returned a non-dict")
    if py_errors:
        report["ok"] = False
        report["errors"] = list(report.get("errors") or []) + py_errors
    if py_warnings:
        report["warnings"] = list(report.get("warnings") or []) + py_warnings
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
