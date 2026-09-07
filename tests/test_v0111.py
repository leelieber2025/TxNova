from __future__ import annotations

import os
from pathlib import Path

import pandas as pd
import pytest
import yaml
from test_pipeline_mock import FIXTURES, FakeAssembler
from txnova.config import QuantifyConfig, load_config, packaged_hexamer_table
from txnova.orchestrator import run_pipeline
from txnova.samples import SampleRow
from txnova.stamps import bam_fingerprint, path_fingerprint
from txnova.stats import run_deseq


def test_fingerprints_detect_same_second_same_size_rewrite(tmp_path):
    path = tmp_path / "sample.bam"
    path.write_bytes(b"old")
    row = SampleRow(sample_id="s", bam=path, strandedness="rf", replicate=1)
    os.utime(path, ns=(1_000_000_001, 1_000_000_001))
    before = path_fingerprint([path]), bam_fingerprint([row])
    path.write_bytes(b"new")
    os.utime(path, ns=(1_000_000_002, 1_000_000_002))
    after = path_fingerprint([path]), bam_fingerprint([row])
    assert before[0] != after[0]
    assert before[1] != after[1]


@pytest.mark.parametrize("change", ["annotation", "fasta", "hexamer", "require_orf", "peptides"])
def test_pipeline_invalidates_changed_inputs(tmp_path, change):
    raw = yaml.safe_load((FIXTURES / "config_ok.yaml").read_text())
    for name in ("genome.fa", "genome.fa.fai", "genes.gtf"):
        (tmp_path / name).write_bytes((FIXTURES / name).read_bytes())
    hexamer = tmp_path / "hexamer.tsv"
    hexamer.write_bytes(packaged_hexamer_table().read_bytes())
    raw.update(output_dir=str(tmp_path / "out"), samples=str(FIXTURES / "samples_ok.tsv"))
    raw["genome"].update(fasta=str(tmp_path / "genome.fa"), annotation=str(tmp_path / "genes.gtf"))
    raw["de"] = {"enabled": False}
    raw["coding"] = {"enabled": True, "fold": False, "orphan": False, "hexamer_table": str(hexamer)}
    cfg_path = tmp_path / "config.yaml"
    cfg_path.write_text(yaml.safe_dump(raw))
    cfg = load_config(cfg_path)
    run_pipeline(cfg, config_path=cfg_path, assembler=FakeAssembler())
    orfs = cfg.output_dir / "candidates/orfs.tsv"
    before = orfs.stat().st_mtime_ns
    if change in {"annotation", "fasta", "hexamer"}:
        path = {"annotation": cfg.genome.annotation, "fasta": cfg.genome.fasta, "hexamer": hexamer}[
            change
        ]
        st = path.stat()
        os.utime(path, ns=(st.st_atime_ns, st.st_mtime_ns + 1_000_000))
    elif change == "require_orf":
        cfg.coding.require_orf = True
    else:
        (cfg.output_dir / "candidates/peptides.fa").unlink()
    run_pipeline(cfg, config_path=cfg_path, assembler=FakeAssembler())
    assert orfs.stat().st_mtime_ns != before
    assert (cfg.output_dir / "candidates/peptides.fa").is_file()
    if change == "annotation":
        import json

        stamp = json.loads((cfg.output_dir / "stamps/assemble.json").read_text())
        assert stamp["fingerprint"]["annotation"] == path_fingerprint([cfg.genome.annotation])


def test_deseq_receives_threshold_and_cpu_limit(tmp_path, monkeypatch):
    from txnova.config import TxNovaConfig

    cfg = TxNovaConfig.model_validate(
        {
            "genome": {"fasta": "/a", "annotation": "/b"},
            "samples": "/s",
            "threads": 2,
            "de": {"padj": 0.01},
        }
    )
    seen = {}

    class DataSet:
        def __init__(self, **kwargs):
            seen["dds"] = kwargs

        def deseq2(self):
            pass

    class Stats:
        def __init__(self, dds, **kwargs):
            seen["stats"] = kwargs
            self.results_df = pd.DataFrame({"padj": [0.001]}, index=["L1"])

        def summary(self):
            pass

    monkeypatch.setattr("pydeseq2.dds.DeseqDataSet", DataSet)
    monkeypatch.setattr("pydeseq2.ds.DeseqStats", Stats)
    rows = [
        SampleRow(sample_id=s, bam=Path("/bam"), group=g, strandedness="rf", replicate=1)
        for s, g in [("c1", "control"), ("c2", "control"), ("t1", "treat"), ("t2", "treat")]
    ]
    counts = tmp_path / "counts.tsv"
    counts.write_text("locus_id\tc1\tc2\tt1\tt2\nL1\t1\t2\t3\t4\n")
    run_deseq(cfg, rows, counts, tmp_path / "de.tsv")
    assert seen["dds"]["n_cpus"] == 2
    assert seen["stats"]["n_cpus"] == 2
    assert seen["stats"]["alpha"] == 0.01


@pytest.mark.parametrize("value", [-1, 256])
def test_invalid_mapq_rejected(value):
    from pydantic import ValidationError

    with pytest.raises(ValidationError):
        QuantifyConfig(min_mapq=value)
