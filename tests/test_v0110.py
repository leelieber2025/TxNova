"""Regression tests for the 0.1.10 review fixes."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
from txnova.config import TxNovaConfig, load_config
from txnova.errors import TxNovaError
from txnova.gates import apply_gates, pick_representative_transcript
from txnova.leak import annotate_leak, exclude_shared_from_finals
from txnova.residual import cluster_junctions
from txnova.rmsk import load_rmsk, rmsk_frac
from txnova.samples import SampleRow
from txnova.stamps import path_fingerprint, stamp_matches, write_stamp

FIXTURES = Path(__file__).resolve().parent / "fixtures"


def _rows_2v2() -> list[SampleRow]:
    return [
        SampleRow(
            sample_id="c1", bam=Path("/c1.bam"), group="control", strandedness="rf", replicate=1
        ),
        SampleRow(
            sample_id="c2", bam=Path("/c2.bam"), group="control", strandedness="rf", replicate=2
        ),
        SampleRow(
            sample_id="t1", bam=Path("/t1.bam"), group="treat", strandedness="rf", replicate=1
        ),
        SampleRow(
            sample_id="t2", bam=Path("/t2.bam"), group="treat", strandedness="rf", replicate=2
        ),
    ]


def _gate_tables(tmp_path: Path, *, control_tpm: float = 0.0, t1: float = 5.0, t2: float = 6.0):
    class_tsv = tmp_path / "c.tsv"
    class_tsv.write_text(
        "transcript_id\tgene_id\tchrom\tstart\tend\tstrand\tn_exons\tlength_nt\texon_structure\tclass\tgene_name\n"
        "T1\tL1\tchr1\t1000\t2000\t+\t2\t400\t1000-1200,1800-2000\tu\tx\n",
        encoding="utf-8",
    )
    reps = tmp_path / "reps.tsv"
    pick_representative_transcript(class_tsv, reps)
    struct = tmp_path / "s.tsv"
    struct.write_text(
        "locus_id\tsample_id\ttranscript_id\tchrom\tstart\tend\tstrand\tn_exons\tlength_nt\t"
        "donors\tacceptors\tcanonical_splice_fraction\tn_introns\thas_nearest\t"
        "nearest_gene_id\tnearest_gene_name\tnearest_distance_bp\tnearest_strand\t"
        "nearest_any_gene_id\tnearest_any_gene_name\tnearest_any_distance_bp\tnearest_any_strand\t"
        "locus_mean_depth\tvalley_mean\tvalley_possible\tgap_mean_depth\tn_dup_flag_seen\t"
        "junction_support\tbridge_read_count\tstructure_error\n"
        "L1\tc1\tT1\tchr1\t1000\t2000\t+\t2\t400\tGT\tAG\t1.0\t1\tfalse\t"
        "\t\t\t\t\t\t\t\t10\t\t\t\t0\t12\t0\n"
        "L1\tt1\tT1\tchr1\t1000\t2000\t+\t2\t400\tGT\tAG\t1.0\t1\tfalse\t"
        "\t\t\t\t\t\t\t\t10\t\t\t\t0\t12\t0\n"
        "L1\tt2\tT1\tchr1\t1000\t2000\t+\t2\t400\tGT\tAG\t1.0\t1\tfalse\t"
        "\t\t\t\t\t\t\t\t10\t\t\t\t0\t8\t0\n",
        encoding="utf-8",
    )
    tpm = tmp_path / "tpm.tsv"
    tpm.write_text(f"locus_id\tc1\tc2\tt1\tt2\nL1\t{control_tpm}\t0.0\t{t1}\t{t2}\n")
    counts = tmp_path / "cnt.tsv"
    counts.write_text("locus_id\tc1\tc2\tt1\tt2\nL1\t0\t0\t10\t12\n")
    return class_tsv, reps, struct, tpm, counts


def test_init_defaults_two_treat_can_pass_detection(tmp_path: Path) -> None:
    cfg = TxNovaConfig.model_validate(
        {"genome": {"fasta": "/a.fa", "annotation": "/a.gtf"}, "samples": "/s.tsv"}
    )
    assert cfg.filters.treat_min_detected_replicates == 2
    class_tsv, reps, struct, tpm, counts = _gate_tables(tmp_path)
    cand = apply_gates(
        cfg,
        _rows_2v2(),
        class_tsv=class_tsv,
        reps_tsv=reps,
        structure_tsv=struct,
        locus_tpm=tpm,
        locus_counts=counts,
        out_dir=tmp_path / "cand",
    )
    assert len(cand) == 1
    assert cand.iloc[0]["locus_id"] == "L1"


def test_rmsk_frac_chrom_key(tmp_path: Path) -> None:
    bed = tmp_path / "rmsk.bed"
    bed.write_text("chr1\t1000\t1100\tAlu\tAlu\n")
    idx = load_rmsk(bed)
    assert abs(rmsk_frac("1", "1001-1200", 200, idx) - 0.5) < 1e-9
    assert abs(rmsk_frac("chr1", "1001-1200", 200, idx) - 0.5) < 1e-9


def test_rmsk_bed_resolves_relative_to_config(tmp_path: Path, monkeypatch) -> None:
    cfg_dir = tmp_path / "cfg"
    other = tmp_path / "cwd"
    cfg_dir.mkdir()
    other.mkdir()
    bed = cfg_dir / "rmsk.bed"
    bed.write_text("chr1\t0\t10\tAlu\tAlu\n")
    sheet = cfg_dir / "s.tsv"
    sheet.write_text("sample_id\tbam\tgroup\tstrandedness\n", encoding="utf-8")
    cfg_path = cfg_dir / "c.yaml"
    cfg_path.write_text(
        "genome:\n  fasta: /a.fa\n  annotation: /a.gtf\n  rmsk_bed: rmsk.bed\nsamples: s.tsv\n",
        encoding="utf-8",
    )
    monkeypatch.chdir(other)
    cfg = load_config(cfg_path)
    assert cfg.genome.rmsk_bed == bed.resolve()


def test_missing_rmsk_bed_is_txnova_error(tmp_path: Path) -> None:
    sheet = tmp_path / "s.tsv"
    sheet.write_text("sample_id\tbam\tgroup\tstrandedness\n", encoding="utf-8")
    p = tmp_path / "c.yaml"
    p.write_text(
        "genome:\n  fasta: /a.fa\n  annotation: /a.gtf\n  rmsk_bed: missing.bed\nsamples: s.tsv\n",
        encoding="utf-8",
    )
    with pytest.raises(TxNovaError, match="rmsk_bed not found"):
        load_config(p)


def test_high_rmsk_control_expressed_absent_from_unnamed(tmp_path: Path) -> None:
    bed = tmp_path / "rmsk.bed"
    bed.write_text("chr1\t999\t2000\tAlu\tAlu\n")
    class_tsv, reps, struct, tpm, counts = _gate_tables(tmp_path, control_tpm=5.0)
    cfg = TxNovaConfig.model_validate(
        {
            "genome": {"fasta": "/a.fa", "annotation": "/a.gtf", "rmsk_bed": str(bed)},
            "samples": "/s.tsv",
        }
    )
    apply_gates(
        cfg,
        _rows_2v2(),
        class_tsv=class_tsv,
        reps_tsv=reps,
        structure_tsv=struct,
        locus_tpm=tpm,
        locus_counts=counts,
        out_dir=tmp_path / "cand",
    )
    unnamed = pd.read_csv(tmp_path / "cand" / "candidates.unnamed.tsv", sep="\t")
    finals = pd.read_csv(tmp_path / "cand" / "candidates.tsv", sep="\t")
    assert unnamed.empty or "L1" not in set(unnamed["locus_id"].astype(str))
    assert finals.empty or "L1" not in set(finals["locus_id"].astype(str))


def test_cluster_plus_not_with_minus() -> None:
    rows = [
        {"chrom": "chr1", "start": 100, "end": 200, "strand": "+"},
        {"chrom": "chr1", "start": 200, "end": 300, "strand": "-"},
    ]
    clusters = cluster_junctions(rows)
    assert len(clusters) == 2


def test_annotate_leak_joins_residual_intron(tmp_path: Path) -> None:
    leak = tmp_path / "leak.tsv"
    leak.write_text(
        "chrom\tstart\tend\tstrand\tstatus\tmerged_locus\nchr1\t451\t499\t+\tunassembled\t\n",
        encoding="utf-8",
    )
    residual = tmp_path / "residual.tsv"
    residual.write_text(
        "residual_id\tchrom\tstrand\tintron_structure\nRSDL.1\tchr1\t+\t451-499\n",
        encoding="utf-8",
    )
    cand = tmp_path / "candidates.tsv"
    cand.write_text("locus_id\nRSDL.1\n", encoding="utf-8")
    gates = tmp_path / "gates.tsv"
    gates.write_text("locus_id\nRSDL.1\n", encoding="utf-8")
    unnamed = tmp_path / "unnamed.tsv"
    unnamed.write_text("locus_id\n", encoding="utf-8")
    df = annotate_leak(leak, gates=gates, unnamed=unnamed, candidates=cand, residual=residual)
    assert bool(df.iloc[0]["in_candidates"]) is True
    assert bool(df.iloc[0]["in_gates"]) is True


def test_gates_stamp_busts_after_bed_change(tmp_path: Path) -> None:
    bed = tmp_path / "rmsk.bed"
    bed.write_text("chr1\t0\t10\n")
    fp1 = {
        "rmsk_bed": str(bed),
        "rmsk": path_fingerprint([bed]),
    }
    stamp = tmp_path / "gates.json"
    write_stamp(stamp, fp1)
    assert stamp_matches(stamp, fp1)
    bed.write_text("chr1\t0\t9999\n")
    fp2 = {"rmsk_bed": str(bed), "rmsk": path_fingerprint([bed])}
    assert not stamp_matches(stamp, fp2)


def test_residual_fp_includes_fai(tmp_path: Path) -> None:
    fai = tmp_path / "g.fa.fai"
    fai.write_text("chr1\t1000\t6\t1000\t1001\n")
    fp1 = {"inputs": path_fingerprint([fai])}
    stamp = tmp_path / "residual.json"
    write_stamp(stamp, fp1)
    fai.write_text("chr1\t2000\t6\t1000\t1001\nchr2\t100\t0\t50\t51\n")
    fp2 = {"inputs": path_fingerprint([fai])}
    assert not stamp_matches(stamp, fp2)


def test_user_agent_is_package_version() -> None:
    from txnova import USER_AGENT, __version__, fold, function, orphan

    assert USER_AGENT == f"txnova/{__version__}"
    assert fold.USER_AGENT == USER_AGENT
    assert "txnova/0.4" not in Path(function.__file__).read_text(encoding="utf-8")
    assert "txnova/0.4" not in Path(orphan.__file__).read_text(encoding="utf-8")


def test_exclude_shared_from_finals_helper(tmp_path: Path) -> None:
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
    n = exclude_shared_from_finals(final, shared, _rows_2v2())
    assert n == 1
    kept = pd.read_csv(final, sep="\t")
    assert list(kept["locus_id"]) == ["RSDL.1"]


def test_orchestrator_calls_exclude_shared() -> None:
    from pathlib import Path as P

    src = P(__file__).resolve().parents[1] / "python" / "txnova" / "orchestrator.py"
    text = src.read_text(encoding="utf-8")
    assert text.count("exclude_shared_from_finals(") >= 3
