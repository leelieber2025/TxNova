from __future__ import annotations

from pathlib import Path

import typer

from txnova import __version__
from txnova.errors import TxNovaError
from txnova.logging import setup_logging

app = typer.Typer(no_args_is_help=True, add_completion=False, pretty_exceptions_enable=False)


def _version_callback(value: bool) -> None:
    if value:
        typer.echo(f"txnova {__version__}")
        raise typer.Exit()


@app.callback()
def _main(
    verbose: bool = typer.Option(False, "--verbose", "-v", help="DEBUG logs"),
    version: bool = typer.Option(
        False,
        "--version",
        callback=_version_callback,
        is_eager=True,
        help="Print version and exit",
    ),
) -> None:
    setup_logging(verbose)


def _fail(exc: Exception) -> None:
    if isinstance(exc, TxNovaError):
        typer.secho(str(exc), err=True, fg=typer.colors.RED)
        raise typer.Exit(code=1)
    raise exc


@app.command()
def init(
    c: Path = typer.Option(..., "-c", "--config", help="Config YAML to write"),
    samples: Path = typer.Option(..., "--samples", help="Sample sheet TSV to write"),
) -> None:
    """Write a starter config.yaml and samples.tsv. Fails if either file exists."""
    from txnova.init_files import write_init_files

    try:
        write_init_files(c, samples)
    except Exception as e:
        _fail(e)
    typer.echo(f"wrote {c} and {samples}")


@app.command()
def preflight(
    c: Path = typer.Option(..., "-c", "--config"),
) -> None:
    """Validate BAMs, FASTA, GTF, and the sample sheet. Writes preflight.json."""
    from txnova.config import load_config
    from txnova.orchestrator import preflight_only

    try:
        cfg = load_config(c)
        preflight_only(cfg)
    except Exception as e:
        _fail(e)
    typer.echo(f"preflight ok → {cfg.output_dir / 'preflight.json'}")


@app.command()
def run(
    c: Path = typer.Option(..., "-c", "--config"),
    force: bool = typer.Option(False, "--force", help="Drop stamps and rerun"),
) -> None:
    """Assemble, classify, quantify, gate, and optionally score DE and coding."""
    from txnova.config import load_config
    from txnova.orchestrator import run_pipeline

    try:
        cfg = load_config(c)
        run_pipeline(cfg, config_path=c, force=force)
    except Exception as e:
        _fail(e)
    cand = cfg.output_dir / "candidates" / "candidates.tsv"
    n = 0
    if cand.is_file():
        n = max(0, sum(1 for _ in cand.open()) - 1)
    typer.echo(f"run complete: {n} candidates under {cfg.output_dir}")


@app.command()
def report(
    c: Path = typer.Option(..., "-c", "--config"),
) -> None:
    """Rewrite Markdown and HTML reports from existing outputs."""
    from txnova.config import load_config
    from txnova.orchestrator import write_report_from_outputs

    try:
        cfg = load_config(c)
        path = write_report_from_outputs(cfg)
    except Exception as e:
        _fail(e)
    typer.echo(str(path))
