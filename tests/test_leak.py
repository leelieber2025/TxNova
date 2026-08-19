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
            "min_samples": 2,
            "min_support": 2,
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
    assert hit.iloc[0]["cohort"] == "silent"


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
            "min_samples": 2,
            "min_support": 2,
            "max_control_support": 0,
        }
    )
    _core.leak_scan(str(merged), payload, str(out), cfg)
    df = pd.read_csv(out, sep="\t")
    hit = df[(df["start"] == 451) & (df["end"] == 499)]
    assert hit.empty


def test_leak_shared_when_control_has_splice(tmp_path: Path) -> None:
    bam = FIXTURES / "junc_bridge.bam"
    merged = tmp_path / "merged.gtf"
    merged.write_text(
        'chr1\tX\ttranscript\t10\t80\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1";\n'
        'chr1\tX\texon\t10\t80\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1";\n',
        encoding="utf-8",
    )
    out = tmp_path / "leak.tsv"
    payload = json.dumps(
        {
            "samples": [
                {"sample_id": "c1", "bam": str(bam), "group": "control", "dup_flag_seen": False},
                {"sample_id": "c2", "bam": str(bam), "group": "control", "dup_flag_seen": False},
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
            "min_samples": 2,
            "min_support": 2,
            "max_control_support": 0,
        }
    )
    _core.leak_scan(str(merged), payload, str(out), cfg)
    df = pd.read_csv(out, sep="\t")
    hit = df[(df["start"] == 451) & (df["end"] == 499)]
    assert len(hit) == 1
    assert hit.iloc[0]["cohort"] == "shared"
    assert int(hit.iloc[0]["control_max"]) >= 1
    assert int(hit.iloc[0]["control_n_detected"]) >= 1
    assert int(hit.iloc[0]["treat_sum"]) >= 2
    assert int(hit.iloc[0]["support_sum"]) >= 2
    assert int(hit.iloc[0]["n_detected"]) >= 2


def test_leak_harvests_without_groups(tmp_path: Path) -> None:
    bam = FIXTURES / "junc_bridge.bam"
    merged = tmp_path / "merged.gtf"
    merged.write_text(
        'chr1\tX\ttranscript\t10\t80\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1";\n'
        'chr1\tX\texon\t10\t80\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1";\n',
        encoding="utf-8",
    )
    out = tmp_path / "leak.tsv"
    payload = json.dumps(
        {
            "samples": [
                {"sample_id": "a", "bam": str(bam), "dup_flag_seen": False},
                {"sample_id": "b", "bam": str(bam), "dup_flag_seen": False},
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
            "min_samples": 2,
            "min_support": 2,
            "max_control_support": 0,
        }
    )
    _core.leak_scan(str(merged), payload, str(out), cfg)
    df = pd.read_csv(out, sep="\t")
    hit = df[(df["start"] == 451) & (df["end"] == 499)]
    assert len(hit) == 1
    assert hit.iloc[0]["cohort"] == "cohort"
    assert int(hit.iloc[0]["n_detected"]) >= 2
    assert int(hit.iloc[0]["support_sum"]) >= 2
    assert int(hit.iloc[0]["treat_sum"]) == 0


def test_write_shared_keeps_control_splice_loci(tmp_path: Path) -> None:
    from txnova.leak import write_shared

    residual = tmp_path / "residual.tsv"
    residual.write_text(
        "residual_id\tcontrol_max\tcontrol_n_detected\tcohort\n"
        "RSDL.1\t0\t0\tsilent\n"
        "RSDL.2\t4\t2\tshared\n",
        encoding="utf-8",
    )
    unnamed = tmp_path / "unnamed.tsv"
    unnamed.write_text(
        "locus_id\tlocus_coord\tn_exons\nRSDL.1\tchr1:1-10:+\t2\nRSDL.2\tchr1:20-40:+\t3\n",
        encoding="utf-8",
    )
    dest = tmp_path / "shared.tsv"
    out = write_shared(residual=residual, structure_tables=[unnamed], dest=dest)
    assert list(out["locus_id"]) == ["RSDL.2"]
    assert str(out.iloc[0]["residual_cohort"]) == "shared"


def test_exclude_shared_from_finals(tmp_path: Path) -> None:
    from txnova.leak import exclude_shared_from_finals
    from txnova.samples import SampleRow

    final = tmp_path / "candidates.tsv"
    final.write_text(
        "locus_id\tlocus_coord\tchrom\tstart\tend\tstrand\tn_exons\tlength_nt\t"
        "exon_structure\tclass\tgates_passed\trepresentative_transcript_id\n"
        "RSDL.1\tchr1:1-10:+\tchr1\t1\t10\t+\t2\t400\t1-10\tu\tclass,detect\tT1\n"
        "RSDL.2\tchr1:20-40:+\tchr1\t20\t40\t+\t3\t400\t20-40\tu\tclass,detect\tT2\n",
        encoding="utf-8",
    )
    shared = tmp_path / "shared.tsv"
    shared.write_text("locus_id\nRSDL.2\n", encoding="utf-8")
    rows = [
        SampleRow(
            sample_id="c1", bam=Path("/c1.bam"), group="control", strandedness="fr", replicate=1
        ),
        SampleRow(
            sample_id="t1", bam=Path("/t1.bam"), group="treat", strandedness="fr", replicate=1
        ),
    ]
    n = exclude_shared_from_finals(final, shared, rows)
    assert n == 1
    kept = [ln.split("\t", 1)[0] for ln in final.read_text().splitlines()[1:] if ln]
    assert kept == ["RSDL.1"]
