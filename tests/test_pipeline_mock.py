from __future__ import annotations

from pathlib import Path

import yaml
from txnova.assemble import Assembler
from txnova.config import load_config
from txnova.orchestrator import run_pipeline
from txnova.stamps import stamp_outputs

FIXTURES = Path(__file__).resolve().parent / "fixtures"

MERGED = """\
chr1\tX\ttranscript\t10\t200\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1";
chr1\tX\texon\t10\t80\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1";
chr1\tX\texon\t120\t200\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1";
chr1\tX\ttranscript\t800\t950\t.\t+\t.\tgene_id "MSTRG.1"; transcript_id "MSTRG.1.1";
chr1\tX\texon\t800\t860\t.\t+\t.\tgene_id "MSTRG.1"; transcript_id "MSTRG.1.1";
chr1\tX\texon\t900\t950\t.\t+\t.\tgene_id "MSTRG.1"; transcript_id "MSTRG.1.1";
"""


class FakeAssembler(Assembler):
    fake_version = "fake-0"

    def resolve_binary(self, name=None):
        return Path("/bin/true")

    def version(self) -> str:
        return self.fake_version

    def assemble(self, bam, gtf_ref, out_gtf, **kwargs) -> None:
        out_gtf.parent.mkdir(parents=True, exist_ok=True)
        out_gtf.write_text(MERGED)
        out_gtf.with_suffix(".log").write_text("fake\n")

    def merge(self, sample_gtfs, gtf_ref, out_gtf, **kwargs) -> None:
        out_gtf.parent.mkdir(parents=True, exist_ok=True)
        out_gtf.write_text(MERGED)
        out_gtf.with_suffix(".log").write_text("fake merge\n")


def test_run_with_mock_assembler(tmp_path: Path) -> None:
    raw = yaml.safe_load((FIXTURES / "config_ok.yaml").read_text())
    raw["output_dir"] = str(tmp_path / "out")
    raw["genome"]["fasta"] = str(FIXTURES / "genome.fa")
    raw["genome"]["annotation"] = str(FIXTURES / "genes.gtf")
    raw["samples"] = str(FIXTURES / "samples_ok.tsv")
    raw["de"] = {"enabled": False}
    raw["coding"] = {"enabled": False}
    raw["filters"] = {
        "class": "u",
        "require_canonical_splice": False,
        "require_coverage_discontinuity": False,
        "treat_median_tpm": 0.0,
        "treat_min_detected_replicates": 1,
        "control_max_tpm": 100.0,
        "transcript_min_nt": 1,
        "min_exons": 2,
    }
    p = tmp_path / "c.yaml"
    p.write_text(yaml.safe_dump(raw))
    cfg = load_config(p)
    run_pipeline(cfg, config_path=p, assembler=FakeAssembler())
    assert (cfg.output_dir / "assembly" / "merged.gtf").is_file()
    assert (cfg.output_dir / "classify" / "transcripts.class.tsv").is_file()
    assert (cfg.output_dir / "quantify" / "locus_tpm.tsv").is_file()
    assert (cfg.output_dir / "candidates" / "candidates.tsv").is_file()
    assert (cfg.output_dir / "candidates" / "residual.tsv").is_file()
    assert (cfg.output_dir / "stamps" / "residual.json").is_file()
    md = (cfg.output_dir / "report" / "report.md").read_text()
    assert "locus/loci in the final table" in md
    assert (cfg.output_dir / "report" / "report.html").is_file()
    leftover = list(cfg.output_dir.glob(".txnova_staging_*"))
    assert leftover == []
    import pandas as pd

    final = pd.read_csv(cfg.output_dir / "candidates" / "candidates.tsv", sep="\t")
    shared_p = cfg.output_dir / "candidates" / "candidates.shared.tsv"
    if shared_p.is_file() and not final.empty and "locus_id" in final.columns:
        shared = pd.read_csv(shared_p, sep="\t")
        if not shared.empty and "locus_id" in shared.columns:
            assert set(final["locus_id"].astype(str)).isdisjoint(
                set(shared["locus_id"].astype(str))
            )
    src = Path(__file__).resolve().parents[1] / "python" / "txnova" / "orchestrator.py"
    assert "exclude_shared_from_finals(" in src.read_text(encoding="utf-8")


def test_second_run_same_outdir_is_idempotent(tmp_path: Path) -> None:
    raw = yaml.safe_load((FIXTURES / "config_ok.yaml").read_text())
    raw["output_dir"] = str(tmp_path / "out")
    raw["genome"]["fasta"] = str(FIXTURES / "genome.fa")
    raw["genome"]["annotation"] = str(FIXTURES / "genes.gtf")
    raw["samples"] = str(FIXTURES / "samples_ok.tsv")
    raw["de"] = {"enabled": False}
    raw["coding"] = {"enabled": False}
    raw["filters"] = {
        "class": "u",
        "require_canonical_splice": False,
        "require_coverage_discontinuity": False,
        "treat_median_tpm": 0.0,
        "treat_min_detected_replicates": 1,
        "control_max_tpm": 100.0,
        "transcript_min_nt": 1,
        "min_exons": 2,
    }
    p = tmp_path / "c.yaml"
    p.write_text(yaml.safe_dump(raw))
    cfg = load_config(p)
    run_pipeline(cfg, config_path=p, assembler=FakeAssembler())
    first = (cfg.output_dir / "candidates" / "candidates.tsv").read_text()
    run_json = (cfg.output_dir / "run.json").read_text()
    assert '"n_merged_transcripts": null' not in run_json
    q_drop = stamp_outputs(cfg.output_dir / "stamps" / "quantify.json").get("dropped_fragments")
    assert q_drop == 0
    run_pipeline(cfg, config_path=p, assembler=FakeAssembler())
    second = (cfg.output_dir / "candidates" / "candidates.tsv").read_text()
    assert first == second
    assert (
        stamp_outputs(cfg.output_dir / "stamps" / "quantify.json").get("dropped_fragments")
        == q_drop
    )
    assert '"dropped_fragments": 0' in (cfg.output_dir / "run.json").read_text()


def test_second_run_does_not_rewrite_orfs(tmp_path: Path) -> None:
    raw = yaml.safe_load((FIXTURES / "config_ok.yaml").read_text())
    raw["output_dir"] = str(tmp_path / "out")
    raw["genome"]["fasta"] = str(FIXTURES / "genome.fa")
    raw["genome"]["annotation"] = str(FIXTURES / "genes.gtf")
    raw["samples"] = str(FIXTURES / "samples_ok.tsv")
    raw["de"] = {"enabled": False}
    raw["coding"] = {"enabled": True, "fold": False, "orphan": False}
    raw["filters"] = {
        "class": "u",
        "require_canonical_splice": False,
        "require_coverage_discontinuity": False,
        "treat_median_tpm": 0.0,
        "treat_min_detected_replicates": 1,
        "control_max_tpm": 100.0,
        "transcript_min_nt": 1,
        "min_exons": 2,
    }
    p = tmp_path / "c.yaml"
    p.write_text(yaml.safe_dump(raw))
    cfg = load_config(p)
    run_pipeline(cfg, config_path=p, assembler=FakeAssembler())
    orfs = cfg.output_dir / "candidates" / "orfs.tsv"
    assert orfs.is_file()
    stamp = (cfg.output_dir / "stamps" / "coding.json").read_text()
    mtime = orfs.stat().st_mtime
    run_pipeline(cfg, config_path=p, assembler=FakeAssembler())
    assert orfs.stat().st_mtime == mtime
    assert (cfg.output_dir / "stamps" / "coding.json").read_text() == stamp
