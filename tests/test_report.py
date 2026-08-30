from __future__ import annotations

from pathlib import Path

import yaml
from txnova.cli import app
from txnova.config import load_config
from txnova.report import render_report
from txnova.samples import load_samples
from typer.testing import CliRunner

FIXTURES = Path(__file__).resolve().parent / "fixtures"


def test_preflight_cli_and_honest_report_text(tmp_path: Path) -> None:
    raw = yaml.safe_load((FIXTURES / "config_ok.yaml").read_text())
    raw["output_dir"] = str(tmp_path / "out")
    raw["genome"]["fasta"] = str(FIXTURES / "genome.fa")
    raw["genome"]["annotation"] = str(FIXTURES / "genes.gtf")
    raw["samples"] = str(FIXTURES / "samples_ok.tsv")
    p = tmp_path / "c.yaml"
    p.write_text(yaml.safe_dump(raw), encoding="utf-8")
    r = CliRunner().invoke(app, ["preflight", "-c", str(p)])
    assert r.exit_code == 0, r.stdout + r.stderr
    cfg = load_config(p)
    rows = load_samples(cfg.samples)
    md = render_report(cfg, rows, {"ok": True, "aligner_family": "STAR"}, n_candidates=0)
    assert "locus/loci in the final table" in md
    assert "0" in md
    assert "## Residual splice loci" in md
    assert "No residual splice loci" in md
    assert "## Gene-like rank" in md
    leftover = list(cfg.output_dir.glob(".txnova_staging_*"))
    assert leftover == []


def test_write_html_report(tmp_path: Path) -> None:
    from txnova.report import render_html, write_report

    raw = yaml.safe_load((FIXTURES / "config_ok.yaml").read_text())
    raw["output_dir"] = str(tmp_path / "out")
    raw["genome"]["fasta"] = str(FIXTURES / "genome.fa")
    raw["genome"]["annotation"] = str(FIXTURES / "genes.gtf")
    raw["samples"] = str(FIXTURES / "samples_ok.tsv")
    p = tmp_path / "c.yaml"
    p.write_text(yaml.safe_dump(raw), encoding="utf-8")
    cfg = load_config(p)
    rows = load_samples(cfg.samples)
    pf = {"ok": True, "aligner_family": "STAR", "warnings": []}
    md = render_report(cfg, rows, pf, n_candidates=0)
    html = render_html(cfg, rows, pf, n_candidates=0)
    write_report(cfg, md, html)
    html_path = cfg.output_dir / "report" / "report.html"
    assert html_path.is_file()
    text = html_path.read_text(encoding="utf-8")
    assert "locus/loci in the final table" in text
    assert "hexamer LLR" in text
    leftover = list(cfg.output_dir.glob(".txnova_staging_*"))
    assert leftover == []
