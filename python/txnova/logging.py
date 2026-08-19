from __future__ import annotations

import logging
import os
import sys


def setup_logging(verbose: bool = False) -> None:
    env = os.environ.get("TXNOVA_LOG", "").strip().lower()
    if env in {"debug", "info", "warning", "error"}:
        level_name = env
    elif verbose:
        level_name = "debug"
    else:
        level_name = "info"
    level = getattr(logging, level_name.upper())
    logging.basicConfig(
        level=level,
        format="%(levelname)s %(name)s: %(message)s",
        stream=sys.stderr,
        force=True,
    )


def get_logger(name: str = "txnova") -> logging.Logger:
    return logging.getLogger(name)
