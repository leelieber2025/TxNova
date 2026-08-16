from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from txnova import _core

FIXTURES = Path(__file__).resolve().parent / "fixtures"


def test_leak_unassembled_treat_specific(tmp_path: Path) -> None:
    bam = FIXTURES / "junc_bridge.bam"
    ctrl = FIXTURES / "ctrl_1.bam"
    merged = tmp_path / "merged.gtf"
    # No intron at 451-499, so the BAM splice is unassembled.
    merged.write_text(
        'chr1\tX\ttranscript\t10\t80\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1";\n'
        'chr1\tX\texon\t10\t80\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1";\n',
        encoding="utf-8",
    )
    out = tmp_path / "leak.tsv"
    payload = json.dumps(
        {
            "samples": [
                {"sample_id": "c1", "bam": str(ctrl), "group": "control", "dup_flag_seen": False},
                {"sample_id": "c2", "bam": str(ctrl), "group": "control", "dup_flag_seen": False},
                {"sample_id": "t1", "bam": str(bam), "group": "treat", "dup_flag_seen": False},
                {"sample_id": "t2", "bam": str(bam), "group": "treat", "dup_flag_seen": False},
            ]
        }
    )
    cfg = json.dumps(
        {
            "strandedness": "unstranded",
            "min_mapq": 10,
            "skip_duplicate": "never",
            "library_layout": "single",
            "threads": 1,
            "treat_min_samples": 2,
            "min_treat_support": 2,
            "max_control_support": 0,
        }
    )
    _core.leak_scan(str(merged), payload, str(out), cfg)
    df = pd.read_csv(out, sep="\t")
    assert not df.empty
    hit = df[(df["start"] == 451) & (df["end"] == 499)]
    assert len(hit) == 1
    assert hit.iloc[0]["status"] == "unassembled"
    assert int(hit.iloc[0]["treat_sum"]) >= 2
    assert int(hit.iloc[0]["control_max"]) == 0


def test_leak_known_intron_omitted(tmp_path: Path) -> None:
    bam = FIXTURES / "junc_bridge.bam"
    ctrl = FIXTURES / "ctrl_1.bam"
    merged = tmp_path / "merged.gtf"
    merged.write_text(
        'chr1\tX\ttranscript\t400\t550\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1";\n'
        'chr1\tX\texon\t400\t450\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1";\n'
        'chr1\tX\texon\t500\t550\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1";\n',
        encoding="utf-8",
    )
    out = tmp_path / "leak.tsv"
    payload = json.dumps(
        {
            "samples": [
                {"sample_id": "c1", "bam": str(ctrl), "group": "control", "dup_flag_seen": False},
                {"sample_id": "t1", "bam": str(bam), "group": "treat", "dup_flag_seen": False},
                {"sample_id": "t2", "bam": str(bam), "group": "treat", "dup_flag_seen": False},
            ]
        }
    )
    cfg = json.dumps(
        {
            "strandedness": "unstranded",
            "min_mapq": 10,
            "skip_duplicate": "never",
            "library_layout": "single",
            "threads": 1,
            "treat_min_samples": 2,
            "min_treat_support": 2,
            "max_control_support": 0,
        }
    )
    _core.leak_scan(str(merged), payload, str(out), cfg)
    df = pd.read_csv(out, sep="\t")
    hit = df[(df["start"] == 451) & (df["end"] == 499)]
    assert hit.empty
