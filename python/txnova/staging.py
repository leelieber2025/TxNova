"""Atomic publish: write under `.txnova_staging_{pid}`, then rename.

Same rules as scdm:
- dead-PID staging dirs are removed on the next run
- a live owner's staging dir is left alone
- failure deletes this process's staging (no half-written finals)
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path

from txnova.logging import get_logger

STAGING_PREFIX = ".txnova_staging_"
log = get_logger("txnova.staging")


def process_is_alive(pid: int) -> bool:
    if pid <= 0:
        return True
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def remove_orphan_staging(outdir: Path) -> list[str]:
    removed: list[str] = []
    if not outdir.is_dir():
        return removed
    for child in outdir.iterdir():
        if not child.is_dir() or not child.name.startswith(STAGING_PREFIX):
            continue
        suffix = child.name[len(STAGING_PREFIX) :]
        try:
            pid = int(suffix)
        except ValueError:
            log.warning("skipping unrecognized staging directory %s", child)
            continue
        if process_is_alive(pid):
            log.warning("keeping staging %s owned by running pid %s", child, pid)
            continue
        shutil.rmtree(child, ignore_errors=True)
        log.warning("removed orphan staging %s (pid %s no longer alive)", child, pid)
        removed.append(child.name)
    return removed


class StagingDir:
    def __init__(self, outdir: Path) -> None:
        self.outdir = outdir
        self.path = outdir / f"{STAGING_PREFIX}{os.getpid()}"

    def __enter__(self) -> StagingDir:
        self.outdir.mkdir(parents=True, exist_ok=True)
        remove_orphan_staging(self.outdir)
        if self.path.exists():
            shutil.rmtree(self.path)
        self.path.mkdir(parents=True)
        return self

    def write_text(self, rel: str, text: str) -> Path:
        dest = self.path / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(text, encoding="utf-8")
        return dest

    def publish(self, names: list[str]) -> None:
        for name in names:
            src = self.path / name
            dst = self.outdir / name
            dst.parent.mkdir(parents=True, exist_ok=True)
            os.replace(src, dst)

    def __exit__(self, exc_type, exc, tb) -> None:
        shutil.rmtree(self.path, ignore_errors=True)
