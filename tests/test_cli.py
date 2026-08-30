from __future__ import annotations

from txnova.cli import app
from typer.testing import CliRunner

runner = CliRunner()


def test_help() -> None:
    r = runner.invoke(app, ["--help"])
    assert r.exit_code == 0
    assert "preflight" in r.stdout
    assert "run" in r.stdout


def test_version() -> None:
    r = runner.invoke(app, ["--version"])
    assert r.exit_code == 0
    assert "txnova" in r.stdout


def test_verbose_preflight_missing_config() -> None:
    r = runner.invoke(app, ["-v", "preflight", "-c", "no_such.yaml"])
    assert r.exit_code == 1
