from __future__ import annotations

from pathlib import Path

import yaml

from txnova.config import load_config
from txnova.preflight import run_preflight, write_preflight_json
from txnova.samples import load_samples

FIXTURES = Path(__file__).resolve().parent / "fixtures"


def _cfg(tmp_path: Path, *, annotation: Path | None = None, sheet: Path | None = None):
    src = FIXTURES / "config_ok.yaml"
    raw = yaml.safe_load(src.read_text())
    raw["output_dir"] = str(tmp_path / "out")
    raw["genome"]["fasta"] = str(FIXTURES / "genome.fa")
    raw["genome"]["annotation"] = str(annotation or (FIXTURES / "genes.gtf"))
    raw["samples"] = str(sheet or (FIXTURES / "samples_ok.tsv"))
    p = tmp_path / "c.yaml"
    p.write_text(yaml.safe_dump(raw), encoding="utf-8")
    return load_config(p)


def test_ok_preflight(tmp_path: Path) -> None:
    cfg = _cfg(tmp_path)
    rows = load_samples(cfg.samples)
    report = run_preflight(cfg, rows)
    assert report["ok"] is True
    assert report["aligner_family"] == "STAR"
    assert report["library_layout"] == "single"
    assert report["n_control"] == 2
    assert report["n_treat"] == 2
    write_preflight_json(cfg, report)
    assert (cfg.output_dir / "preflight.json").is_file()
    leftover = list(cfg.output_dir.glob(".txnova_staging_*"))
    assert leftover == []


def test_fail_gtf_chr_prefix(tmp_path: Path) -> None:
    cfg = _cfg(tmp_path, annotation=FIXTURES / "genes_nochr.gtf")
    rows = load_samples(cfg.samples)
    report = run_preflight(cfg, rows)
    assert report["ok"] is False
    assert any("chr" in e.lower() for e in report["errors"])


def test_fail_missing_fai(tmp_path: Path) -> None:
    fa = tmp_path / "genome.fa"
    fa.write_text((FIXTURES / "genome.fa").read_text())
    cfg = _cfg(tmp_path)
    cfg.genome.fasta = fa
    rows = load_samples(cfg.samples)
    report = run_preflight(cfg, rows)
    assert report["ok"] is False
    assert any("faidx" in e or ".fai" in e for e in report["errors"])


def test_fail_truncated(tmp_path: Path) -> None:
    sheet = tmp_path / "s.tsv"
    sheet.write_text(
        "sample_id\tbam\tgroup\tstrandedness\treplicate\n"
        f"ctrl_1\t{FIXTURES / 'truncated.bam'}\tcontrol\trf\t1\n"
        f"ctrl_2\t{FIXTURES / 'ctrl_2.bam'}\tcontrol\trf\t2\n"
        f"treat_1\t{FIXTURES / 'treat_1.bam'}\ttreat\trf\t1\n"
        f"treat_2\t{FIXTURES / 'treat_2.bam'}\ttreat\trf\t2\n",
        encoding="utf-8",
    )
    cfg = _cfg(tmp_path, sheet=sheet)
    rows = load_samples(cfg.samples)
    report = run_preflight(cfg, rows)
    assert report["ok"] is False
    assert any("EOF" in e or "truncated" in e.lower() for e in report["errors"])


def test_fail_no_index(tmp_path: Path) -> None:
    sheet = tmp_path / "s.tsv"
    sheet.write_text(
        "sample_id\tbam\tgroup\tstrandedness\treplicate\n"
        f"ctrl_1\t{FIXTURES / 'noindex.bam'}\tcontrol\trf\t1\n"
        f"ctrl_2\t{FIXTURES / 'ctrl_2.bam'}\tcontrol\trf\t2\n"
        f"treat_1\t{FIXTURES / 'treat_1.bam'}\ttreat\trf\t1\n"
        f"treat_2\t{FIXTURES / 'treat_2.bam'}\ttreat\trf\t2\n",
        encoding="utf-8",
    )
    cfg = _cfg(tmp_path, sheet=sheet)
    rows = load_samples(cfg.samples)
    report = run_preflight(cfg, rows)
    assert report["ok"] is False


def test_fail_bowtie2(tmp_path: Path) -> None:
    sheet = tmp_path / "s.tsv"
    sheet.write_text(
        "sample_id\tbam\tgroup\tstrandedness\treplicate\n"
        f"ctrl_1\t{FIXTURES / 'bowtie2.bam'}\tcontrol\trf\t1\n"
        f"ctrl_2\t{FIXTURES / 'ctrl_2.bam'}\tcontrol\trf\t2\n"
        f"treat_1\t{FIXTURES / 'treat_1.bam'}\ttreat\trf\t1\n"
        f"treat_2\t{FIXTURES / 'treat_2.bam'}\ttreat\trf\t2\n",
        encoding="utf-8",
    )
    cfg = _cfg(tmp_path, sheet=sheet)
    rows = load_samples(cfg.samples)
    report = run_preflight(cfg, rows)
    assert report["ok"] is False
    assert any("Bowtie2" in e for e in report["errors"])


def test_fail_star_plus_bowtie2(tmp_path: Path) -> None:
    sheet = tmp_path / "s.tsv"
    sheet.write_text(
        "sample_id\tbam\tgroup\tstrandedness\treplicate\n"
        f"ctrl_1\t{FIXTURES / 'star_bowtie2.bam'}\tcontrol\trf\t1\n"
        f"ctrl_2\t{FIXTURES / 'ctrl_2.bam'}\tcontrol\trf\t2\n"
        f"treat_1\t{FIXTURES / 'treat_1.bam'}\ttreat\trf\t1\n"
        f"treat_2\t{FIXTURES / 'treat_2.bam'}\ttreat\trf\t2\n",
        encoding="utf-8",
    )
    cfg = _cfg(tmp_path, sheet=sheet)
    rows = load_samples(cfg.samples)
    report = run_preflight(cfg, rows)
    assert report["ok"] is False
    assert any("Bowtie2" in e for e in report["errors"])


def test_fail_mixed_strandedness(tmp_path: Path) -> None:
    sheet = tmp_path / "s.tsv"
    sheet.write_text(
        "sample_id\tbam\tgroup\tstrandedness\treplicate\n"
        f"ctrl_1\t{FIXTURES / 'ctrl_1.bam'}\tcontrol\trf\t1\n"
        f"ctrl_2\t{FIXTURES / 'ctrl_2.bam'}\tcontrol\tfr\t2\n"
        f"treat_1\t{FIXTURES / 'treat_1.bam'}\ttreat\trf\t1\n"
        f"treat_2\t{FIXTURES / 'treat_2.bam'}\ttreat\trf\t2\n",
        encoding="utf-8",
    )
    cfg = _cfg(tmp_path, sheet=sheet)
    rows = load_samples(cfg.samples)
    report = run_preflight(cfg, rows)
    assert report["ok"] is False
    assert any("strandedness" in e for e in report["errors"])


def test_cli_preflight_ok(tmp_path: Path) -> None:
    from typer.testing import CliRunner

    from txnova.cli import app

    cfg = _cfg(tmp_path)
    r = CliRunner().invoke(app, ["preflight", "-c", str(tmp_path / "c.yaml")])
    assert r.exit_code == 0, r.stdout + r.stderr
    assert (cfg.output_dir / "preflight.json").is_file()


def test_ok_without_groups(tmp_path: Path) -> None:
    sheet = tmp_path / "s.tsv"
    sheet.write_text(
        "sample_id\tbam\tstrandedness\n"
        f"a\t{FIXTURES / 'treat_1.bam'}\trf\n"
        f"b\t{FIXTURES / 'treat_2.bam'}\trf\n",
        encoding="utf-8",
    )
    cfg = _cfg(tmp_path, sheet=sheet)
    rows = load_samples(cfg.samples)
    report = run_preflight(cfg, rows)
    assert report["ok"] is True, report.get("errors")
    assert report["n_control"] == 0
    assert report["n_treat"] == 0
    assert cfg.de.enabled is False


def test_ok_treat_only(tmp_path: Path) -> None:
    sheet = tmp_path / "s.tsv"
    sheet.write_text(
        "sample_id\tbam\tgroup\tstrandedness\n"
        f"t1\t{FIXTURES / 'treat_1.bam'}\ttreat\trf\n"
        f"t2\t{FIXTURES / 'treat_2.bam'}\ttreat\trf\n",
        encoding="utf-8",
    )
    cfg = _cfg(tmp_path, sheet=sheet)
    rows = load_samples(cfg.samples)
    report = run_preflight(cfg, rows)
    assert report["ok"] is True, report.get("errors")
    assert report["n_treat"] == 2
    assert report["n_control"] == 0
    assert cfg.de.enabled is False
