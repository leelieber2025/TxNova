from __future__ import annotations

from pathlib import Path

import pandas as pd

from txnova.config import TxNovaConfig
from txnova.gates import (
    STRUCTURE_FEATURE_COLUMNS,
    apply_gates,
    is_structure_error,
    pick_representative_transcript,
    validate_structure_features,
)
from txnova.samples import SampleRow


def test_validate_structure_features_requires_full_schema_and_aligned_introns() -> None:
    from txnova.errors import TxNovaError

    df = pd.DataFrame([{c: "" for c in STRUCTURE_FEATURE_COLUMNS}], dtype=object)
    df.loc[0, "locus_id"] = "L1"
    df.loc[0, "sample_id"] = "t1"
    df.loc[0, "n_introns"] = 1
    df.loc[0, "donors"] = "GT"
    df.loc[0, "acceptors"] = "AG"
    df.loc[0, "junction_support"] = "4"
    validate_structure_features(df, Path("structure.features.tsv"))

    short = df.drop(columns=["structure_error"])
    try:
        validate_structure_features(short, Path("old.tsv"))
    except TxNovaError as e:
        assert "structure_error" in str(e)
    else:
        raise AssertionError("expected missing column")

    bad = df.copy()
    bad.loc[0, "junction_support"] = "4,1"
    try:
        validate_structure_features(bad, Path("misaligned.tsv"))
    except TxNovaError as e:
        assert "junction_support" in str(e)
    else:
        raise AssertionError("expected length mismatch")


def test_is_structure_error() -> None:
    assert not is_structure_error("")
    assert not is_structure_error("nan")
    assert not is_structure_error(None)
    assert not is_structure_error(float("nan"))
    assert is_structure_error("overlapping exons")


def test_pick_longest_then_id(tmp_path: Path) -> None:
    p = tmp_path / "c.tsv"
    p.write_text(
        "transcript_id\tgene_id\tchrom\tstart\tend\tstrand\tn_exons\tlength_nt\texon_structure\tclass\tgene_name\n"
        "M.1.b\tM.1\tchr1\t1\t100\t+\t2\t80\t1-40,60-100\tu\tM\n"
        "M.1.a\tM.1\tchr1\t1\t90\t+\t2\t80\t1-40,50-90\tu\tM\n"
        "M.1.c\tM.1\tchr1\t1\t50\t+\t2\t40\t1-50\tu\tM\n",
        encoding="utf-8",
    )
    reps = pick_representative_transcript(p, tmp_path / "reps.tsv")
    assert list(reps["transcript_id"]) == ["M.1.a"]


def test_detect_gate(tmp_path: Path) -> None:
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
    tpm.write_text("locus_id\tc1\tc2\tt1\tt2\nL1\t0.0\t0.0\t5.0\t6.0\n")
    counts = tmp_path / "cnt.tsv"
    counts.write_text("locus_id\tc1\tc2\tt1\tt2\nL1\t0\t0\t10\t12\n")
    cfg = TxNovaConfig.model_validate(
        {"genome": {"fasta": "/a.fa", "annotation": "/a.gtf"}, "samples": "/s.tsv"}
    )
    rows = [
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
    cand = apply_gates(
        cfg,
        rows,
        class_tsv=class_tsv,
        reps_tsv=reps,
        structure_tsv=struct,
        locus_tpm=tpm,
        locus_counts=counts,
        out_dir=tmp_path / "cand",
    )
    assert len(cand) == 1
    assert cand.iloc[0]["locus_id"] == "L1"
    assert cand.iloc[0]["exon_structure"] == "1000-1200,1800-2000"
    assert cand.iloc[0]["junction_support"] == "20"
    assert int(cand.iloc[0]["junction_support_min"]) == 20
    assert int(cand.iloc[0]["bridge_read_count"]) == 0
    assert "no_bridge" in str(cand.iloc[0]["gates_passed"])
    assert "transcript_min_nt" in str(cand.iloc[0]["gates_passed"])
    assert "min_exons" in str(cand.iloc[0]["gates_passed"])
    assert "multi_exon" not in str(cand.iloc[0]["gates_passed"])
    unnamed = tmp_path / "cand" / "candidates.unnamed.tsv"
    assert unnamed.is_file()
    assert max(0, sum(1 for _ in unnamed.open()) - 1) == 0


def test_control_expressed_goes_to_unnamed(tmp_path: Path) -> None:
    class_tsv = tmp_path / "c.tsv"
    class_tsv.write_text(
        "transcript_id\tgene_id\tchrom\tstart\tend\tstrand\tn_exons\tlength_nt\texon_structure\tclass\tgene_name\n"
        "T1\tL1\tchr1\t1000\t2000\t+\t2\t400\t1000-1200,1800-2000\tu\t\n",
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
    tpm.write_text("locus_id\tc1\tc2\tt1\tt2\nL1\t2.0\t3.0\t5.0\t6.0\n")
    counts = tmp_path / "cnt.tsv"
    counts.write_text("locus_id\tc1\tc2\tt1\tt2\nL1\t4\t6\t10\t12\n")
    cfg = TxNovaConfig.model_validate(
        {"genome": {"fasta": "/a.fa", "annotation": "/a.gtf"}, "samples": "/s.tsv"}
    )
    rows = [
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
    cand = apply_gates(
        cfg,
        rows,
        class_tsv=class_tsv,
        reps_tsv=reps,
        structure_tsv=struct,
        locus_tpm=tpm,
        locus_counts=counts,
        out_dir=tmp_path / "cand",
    )
    assert cand.empty
    unnamed = pd.read_csv(tmp_path / "cand" / "candidates.unnamed.tsv", sep="\t")
    assert len(unnamed) == 1
    assert unnamed.iloc[0]["locus_id"] == "L1"
    assert float(unnamed.iloc[0]["control_max_tpm"]) >= 1.0
    assert "detect" not in str(unnamed.iloc[0]["gates_passed"])


def test_bridge_rejects_readthrough(tmp_path: Path) -> None:
    class_tsv = tmp_path / "c.tsv"
    class_tsv.write_text(
        "transcript_id\tgene_id\tchrom\tstart\tend\tstrand\tn_exons\tlength_nt\texon_structure\tclass\tgene_name\n"
        "T1\tL1\tchr1\t4000\t5000\t+\t2\t400\t4000-4200,4800-5000\tu\tx\n",
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
        "L1\tc1\tT1\tchr1\t4000\t5000\t+\t2\t400\tGT\tAG\t1.0\t1\ttrue\t"
        "G1\tCar\t5416\t+\tG1\tCar\t5416\t+\t10\t0.0\ttrue\t0.0\t0\t0\t0\n"
        "L1\tt1\tT1\tchr1\t4000\t5000\t+\t2\t400\tGT\tAG\t1.0\t1\ttrue\t"
        "G1\tCar\t5416\t+\tG1\tCar\t5416\t+\t10\t0.0\ttrue\t0.0\t0\t20\t3\n"
        "L1\tt2\tT1\tchr1\t4000\t5000\t+\t2\t400\tGT\tAG\t1.0\t1\ttrue\t"
        "G1\tCar\t5416\t+\tG1\tCar\t5416\t+\t10\t0.0\ttrue\t0.0\t0\t18\t2\n",
        encoding="utf-8",
    )
    tpm = tmp_path / "tpm.tsv"
    tpm.write_text("locus_id\tc1\tc2\tt1\tt2\nL1\t0.0\t0.0\t5.0\t6.0\n")
    counts = tmp_path / "cnt.tsv"
    counts.write_text("locus_id\tc1\tc2\tt1\tt2\nL1\t0\t0\t10\t12\n")
    cfg = TxNovaConfig.model_validate(
        {"genome": {"fasta": "/a.fa", "annotation": "/a.gtf"}, "samples": "/s.tsv"}
    )
    rows = [
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
    cand = apply_gates(
        cfg,
        rows,
        class_tsv=class_tsv,
        reps_tsv=reps,
        structure_tsv=struct,
        locus_tpm=tpm,
        locus_counts=counts,
        out_dir=tmp_path / "cand",
    )
    assert cand.empty


def _minimal_gate_inputs(tmp_path: Path, treat_tpms: str, has_nearest: str, dist: str):
    class_tsv = tmp_path / "c.tsv"
    class_tsv.write_text(
        "transcript_id\tgene_id\tchrom\tstart\tend\tstrand\tn_exons\tlength_nt\texon_structure\tclass\tgene_name\n"
        "T1\tL1\tchr1\t1000\t2000\t+\t2\t400\t1000-1200,1800-2000\tu\tx\n",
        encoding="utf-8",
    )
    reps = tmp_path / "reps.tsv"
    pick_representative_transcript(class_tsv, reps)
    struct = tmp_path / "s.tsv"
    nearest_cols = (
        f"{has_nearest}\tG1\tCar\t{dist}\t+\tG1\tCar\t{dist}\t+"
        if has_nearest == "true"
        else "false\t\t\t\t\t\t\t\t"
    )
    struct.write_text(
        "locus_id\tsample_id\ttranscript_id\tchrom\tstart\tend\tstrand\tn_exons\tlength_nt\t"
        "donors\tacceptors\tcanonical_splice_fraction\tn_introns\thas_nearest\t"
        "nearest_gene_id\tnearest_gene_name\tnearest_distance_bp\tnearest_strand\t"
        "nearest_any_gene_id\tnearest_any_gene_name\tnearest_any_distance_bp\tnearest_any_strand\t"
        "locus_mean_depth\tvalley_mean\tvalley_possible\tgap_mean_depth\tn_dup_flag_seen\t"
        "junction_support\tbridge_read_count\tstructure_error\n"
        f"L1\tc1\tT1\tchr1\t1000\t2000\t+\t2\t400\tGT\tAG\t1.0\t1\t{nearest_cols}\t"
        "10\t0.0\ttrue\t0.0\t0\t12\t0\n"
        f"L1\tt1\tT1\tchr1\t1000\t2000\t+\t2\t400\tGT\tAG\t1.0\t1\t{nearest_cols}\t"
        "10\t0.0\ttrue\t0.0\t0\t12\t0\n"
        f"L1\tt2\tT1\tchr1\t1000\t2000\t+\t2\t400\tGT\tAG\t1.0\t1\t{nearest_cols}\t"
        "10\t0.0\ttrue\t0.0\t0\t8\t0\n",
        encoding="utf-8",
    )
    tpm = tmp_path / "tpm.tsv"
    tpm.write_text(f"locus_id\tc1\tc2\tt1\tt2\nL1\t0.0\t0.0\t{treat_tpms}\n")
    counts = tmp_path / "cnt.tsv"
    counts.write_text("locus_id\tc1\tc2\tt1\tt2\nL1\t0\t0\t10\t12\n")
    cfg = TxNovaConfig.model_validate(
        {"genome": {"fasta": "/a.fa", "annotation": "/a.gtf"}, "samples": "/s.tsv"}
    )
    cfg.filters.require_coverage_discontinuity = False
    cfg.filters.min_nearest_same_strand_bp = 0
    cfg.filters.treat_median_tpm = 1.0
    cfg.filters.control_max_tpm = 0.1
    cfg.filters.transcript_min_nt = 200
    rows = [
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
    return cfg, rows, class_tsv, reps, struct, tpm, counts


def test_even_n_median_is_not_upper_half(tmp_path: Path) -> None:
    # treat TPMs 0.5, 1.2 → true median 0.85 < treat_median_tpm 1.0
    # the old [n//2] pick was 1.2 and would have passed
    cfg, rows, class_tsv, reps, struct, tpm, counts = _minimal_gate_inputs(
        tmp_path, "0.5\t1.2", "false", ""
    )
    cand = apply_gates(
        cfg,
        rows,
        class_tsv=class_tsv,
        reps_tsv=reps,
        structure_tsv=struct,
        locus_tpm=tpm,
        locus_counts=counts,
        out_dir=tmp_path / "cand",
    )
    assert cand.empty


def test_has_nearest_without_distance_is_error(tmp_path: Path) -> None:
    from txnova.errors import TxNovaError

    cfg, rows, class_tsv, reps, struct, tpm, counts = _minimal_gate_inputs(
        tmp_path, "5.0\t6.0", "true", ""
    )
    try:
        apply_gates(
            cfg,
            rows,
            class_tsv=class_tsv,
            reps_tsv=reps,
            structure_tsv=struct,
            locus_tpm=tpm,
            locus_counts=counts,
            out_dir=tmp_path / "cand",
        )
    except TxNovaError as e:
        assert "has_nearest=true" in str(e)
    else:
        raise AssertionError("expected hard error for has_nearest without distance")


def test_empty_candidates_have_full_schema(tmp_path: Path) -> None:
    class_tsv = tmp_path / "c.tsv"
    class_tsv.write_text(
        "transcript_id\tgene_id\tchrom\tstart\tend\tstrand\tn_exons\tlength_nt\texon_structure\tclass\tgene_name\n"
        "T1\tL1\tchr1\t1000\t2000\t+\t1\t100\t1000-2000\tu\tx\n",
        encoding="utf-8",
    )
    reps = tmp_path / "reps.tsv"
    pick_representative_transcript(class_tsv, reps)
    cfg = TxNovaConfig.model_validate(
        {"genome": {"fasta": "/a.fa", "annotation": "/a.gtf"}, "samples": "/s.tsv"}
    )
    rows = [
        SampleRow(
            sample_id="c1", bam=Path("/c1.bam"), group="control", strandedness="rf", replicate=1
        ),
        SampleRow(
            sample_id="t1", bam=Path("/t1.bam"), group="treat", strandedness="rf", replicate=1
        ),
        SampleRow(
            sample_id="t2", bam=Path("/t2.bam"), group="treat", strandedness="rf", replicate=2
        ),
    ]
    tpm = tmp_path / "tpm.tsv"
    tpm.write_text("locus_id\tc1\tt1\tt2\nL1\t0.0\t0.0\t0.0\n")
    counts = tmp_path / "cnt.tsv"
    counts.write_text("locus_id\tc1\tt1\tt2\nL1\t0\t0\t0\n")
    cand = apply_gates(
        cfg,
        rows,
        class_tsv=class_tsv,
        reps_tsv=reps,
        structure_tsv=tmp_path / "missing.tsv",
        locus_tpm=tpm,
        locus_counts=counts,
        out_dir=tmp_path / "cand",
    )
    assert cand.empty
    written = pd.read_csv(tmp_path / "cand" / "candidates.tsv", sep="\t")
    assert written.empty
    for col in (
        "locus_id",
        "gates_passed",
        "treat_median_tpm",
        "c1_tpm",
        "t1_count",
        "padj",
        "coding_label",
    ):
        assert col in written.columns


def test_missing_tpm_locus_is_error(tmp_path: Path) -> None:
    from txnova.errors import TxNovaError

    class_tsv = tmp_path / "c.tsv"
    class_tsv.write_text(
        "transcript_id\tgene_id\tchrom\tstart\tend\tstrand\tn_exons\tlength_nt\texon_structure\tclass\tgene_name\n"
        "T1\tL1\tchr1\t1000\t2000\t+\t2\t400\t1000-1200,1800-2000\tu\tx\n",
        encoding="utf-8",
    )
    reps = tmp_path / "reps.tsv"
    pick_representative_transcript(class_tsv, reps)
    tpm = tmp_path / "tpm.tsv"
    tpm.write_text("locus_id\tc1\tt1\nOTHER\t0.0\t1.0\n")
    cfg = TxNovaConfig.model_validate(
        {"genome": {"fasta": "/a.fa", "annotation": "/a.gtf"}, "samples": "/s.tsv"}
    )
    rows = [
        SampleRow(
            sample_id="c1", bam=Path("/c1.bam"), group="control", strandedness="rf", replicate=1
        ),
        SampleRow(
            sample_id="t1", bam=Path("/t1.bam"), group="treat", strandedness="rf", replicate=1
        ),
    ]
    try:
        apply_gates(
            cfg,
            rows,
            class_tsv=class_tsv,
            reps_tsv=reps,
            structure_tsv=tmp_path / "missing.tsv",
            locus_tpm=tpm,
            locus_counts=tmp_path / "cnt.tsv",
            out_dir=tmp_path / "cand",
        )
    except TxNovaError as e:
        assert "missing locus_id L1" in str(e)
    else:
        raise AssertionError("expected missing TPM locus to fail closed")


def test_missing_gap_mean_column_is_error(tmp_path: Path) -> None:
    from txnova.errors import TxNovaError

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
        "locus_mean_depth\tvalley_mean\tvalley_possible\tn_dup_flag_seen\t"
        "junction_support\tbridge_read_count\tstructure_error\n"
        "L1\tt1\tT1\tchr1\t1000\t2000\t+\t2\t400\tGT\tAG\t1.0\t1\tfalse\t"
        "\t\t\t\t\t\t\t\t10\t\t\t0\t12\t0\n",
        encoding="utf-8",
    )
    tpm = tmp_path / "tpm.tsv"
    tpm.write_text("locus_id\tc1\tt1\tt2\nL1\t0.0\t5.0\t6.0\n")
    counts = tmp_path / "cnt.tsv"
    counts.write_text("locus_id\tc1\tt1\tt2\nL1\t0\t10\t12\n")
    cfg = TxNovaConfig.model_validate(
        {"genome": {"fasta": "/a.fa", "annotation": "/a.gtf"}, "samples": "/s.tsv"}
    )
    rows = [
        SampleRow(
            sample_id="c1", bam=Path("/c1.bam"), group="control", strandedness="rf", replicate=1
        ),
        SampleRow(
            sample_id="t1", bam=Path("/t1.bam"), group="treat", strandedness="rf", replicate=1
        ),
        SampleRow(
            sample_id="t2", bam=Path("/t2.bam"), group="treat", strandedness="rf", replicate=2
        ),
    ]
    try:
        apply_gates(
            cfg,
            rows,
            class_tsv=class_tsv,
            reps_tsv=reps,
            structure_tsv=struct,
            locus_tpm=tpm,
            locus_counts=counts,
            out_dir=tmp_path / "cand",
        )
    except TxNovaError as e:
        assert "gap_mean_depth" in str(e)
    else:
        raise AssertionError("expected missing gap_mean_depth to fail closed")


def test_missing_structure_error_column_is_error(tmp_path: Path) -> None:
    from txnova.errors import TxNovaError

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
        "junction_support\tbridge_read_count\n"
        "L1\tt1\tT1\tchr1\t1000\t2000\t+\t2\t400\tGT\tAG\t1.0\t1\tfalse\t"
        "\t\t\t\t\t\t\t\t10\t\t\t0.0\t0\t12\t0\n",
        encoding="utf-8",
    )
    tpm = tmp_path / "tpm.tsv"
    tpm.write_text("locus_id\tc1\tt1\tt2\nL1\t0.0\t5.0\t6.0\n")
    counts = tmp_path / "cnt.tsv"
    counts.write_text("locus_id\tc1\tt1\tt2\nL1\t0\t10\t12\n")
    cfg = TxNovaConfig.model_validate(
        {"genome": {"fasta": "/a.fa", "annotation": "/a.gtf"}, "samples": "/s.tsv"}
    )
    rows = [
        SampleRow(
            sample_id="c1", bam=Path("/c1.bam"), group="control", strandedness="rf", replicate=1
        ),
        SampleRow(
            sample_id="t1", bam=Path("/t1.bam"), group="treat", strandedness="rf", replicate=1
        ),
        SampleRow(
            sample_id="t2", bam=Path("/t2.bam"), group="treat", strandedness="rf", replicate=2
        ),
    ]
    try:
        apply_gates(
            cfg,
            rows,
            class_tsv=class_tsv,
            reps_tsv=reps,
            structure_tsv=struct,
            locus_tpm=tpm,
            locus_counts=counts,
            out_dir=tmp_path / "cand",
        )
    except TxNovaError as e:
        assert "structure_error" in str(e)
    else:
        raise AssertionError("expected missing structure_error to fail closed")


def test_junction_support_length_must_match_n_introns(tmp_path: Path) -> None:
    from txnova.errors import TxNovaError

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
        "\t\t\t\t\t\t\t\t10\t\t\t\t0\t12,3\t0\t\n"
        "L1\tt1\tT1\tchr1\t1000\t2000\t+\t2\t400\tGT\tAG\t1.0\t1\tfalse\t"
        "\t\t\t\t\t\t\t\t10\t\t\t\t0\t12,3\t0\t\n"
        "L1\tt2\tT1\tchr1\t1000\t2000\t+\t2\t400\tGT\tAG\t1.0\t1\tfalse\t"
        "\t\t\t\t\t\t\t\t10\t\t\t\t0\t8,1\t0\t\n",
        encoding="utf-8",
    )
    tpm = tmp_path / "tpm.tsv"
    tpm.write_text("locus_id\tc1\tc2\tt1\tt2\nL1\t0.0\t0.0\t5.0\t6.0\n")
    counts = tmp_path / "cnt.tsv"
    counts.write_text("locus_id\tc1\tc2\tt1\tt2\nL1\t0\t0\t10\t12\n")
    cfg = TxNovaConfig.model_validate(
        {"genome": {"fasta": "/a.fa", "annotation": "/a.gtf"}, "samples": "/s.tsv"}
    )
    rows = [
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
    try:
        apply_gates(
            cfg,
            rows,
            class_tsv=class_tsv,
            reps_tsv=reps,
            structure_tsv=struct,
            locus_tpm=tpm,
            locus_counts=counts,
            out_dir=tmp_path / "cand",
        )
    except TxNovaError as e:
        assert "junction_support" in str(e) and "n_introns" in str(e)
    else:
        raise AssertionError("expected junction_support / n_introns mismatch to fail closed")


def test_transcribed_gap_rejected_despite_empty_valley(tmp_path: Path) -> None:
    """9252-class: 200 bp valley empty, but the gap as a whole is transcribed."""
    class_tsv = tmp_path / "c.tsv"
    class_tsv.write_text(
        "transcript_id\tgene_id\tchrom\tstart\tend\tstrand\tn_exons\tlength_nt\texon_structure\tclass\tgene_name\n"
        "T1\tL1\tchr3\t14733479\t14738935\t+\t2\t4771\t14733479-14734062,14734749-14738935\tu\tx\n",
        encoding="utf-8",
    )
    reps = tmp_path / "reps.tsv"
    pick_representative_transcript(class_tsv, reps)
    header = (
        "locus_id\tsample_id\ttranscript_id\tchrom\tstart\tend\tstrand\tn_exons\tlength_nt\t"
        "donors\tacceptors\tcanonical_splice_fraction\tn_introns\thas_nearest\t"
        "nearest_gene_id\tnearest_gene_name\tnearest_distance_bp\tnearest_strand\t"
        "nearest_any_gene_id\tnearest_any_gene_name\tnearest_any_distance_bp\tnearest_any_strand\t"
        "locus_mean_depth\tvalley_mean\tvalley_possible\tgap_mean_depth\tn_dup_flag_seen\t"
        "junction_support\tbridge_read_count\tstructure_error\n"
    )
    struct = tmp_path / "s.tsv"
    struct.write_text(
        header + "L1\tc1\tT1\tchr3\t14733479\t14738935\t+\t2\t4771\tGT\tAG\t1.0\t1\ttrue\t"
        "G1\tCar13\t5416\t+\tG1\tCar13\t5416\t+\t10\t0.0\ttrue\t12.2\t0\t0\t0\n"
        "L1\tt1\tT1\tchr3\t14733479\t14738935\t+\t2\t4771\tGT\tAG\t1.0\t1\ttrue\t"
        "G1\tCar13\t5416\t+\tG1\tCar13\t5416\t+\t11\t0.0\ttrue\t51.0\t0\t1\t0\n"
        "L1\tt2\tT1\tchr3\t14733479\t14738935\t+\t2\t4771\tGT\tAG\t1.0\t1\ttrue\t"
        "G1\tCar13\t5416\t+\tG1\tCar13\t5416\t+\t20\t0.0\ttrue\t63.8\t0\t1\t0\n",
        encoding="utf-8",
    )
    tpm = tmp_path / "tpm.tsv"
    tpm.write_text("locus_id\tc1\tc2\tt1\tt2\nL1\t0.0\t0.0\t5.0\t6.0\n")
    counts = tmp_path / "cnt.tsv"
    counts.write_text("locus_id\tc1\tc2\tt1\tt2\nL1\t0\t0\t10\t12\n")
    cfg = TxNovaConfig.model_validate(
        {"genome": {"fasta": "/a.fa", "annotation": "/a.gtf"}, "samples": "/s.tsv"}
    )
    rows = [
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
    cand = apply_gates(
        cfg,
        rows,
        class_tsv=class_tsv,
        reps_tsv=reps,
        structure_tsv=struct,
        locus_tpm=tpm,
        locus_counts=counts,
        out_dir=tmp_path / "cand",
    )
    assert cand.empty


def test_empty_intergenic_gap_passes(tmp_path: Path) -> None:
    class_tsv = tmp_path / "c.tsv"
    class_tsv.write_text(
        "transcript_id\tgene_id\tchrom\tstart\tend\tstrand\tn_exons\tlength_nt\texon_structure\tclass\tgene_name\n"
        "T1\tL1\tchr5\t1000\t2000\t+\t2\t400\t1000-1200,1800-2000\tu\tx\n",
        encoding="utf-8",
    )
    reps = tmp_path / "reps.tsv"
    pick_representative_transcript(class_tsv, reps)
    header = (
        "locus_id\tsample_id\ttranscript_id\tchrom\tstart\tend\tstrand\tn_exons\tlength_nt\t"
        "donors\tacceptors\tcanonical_splice_fraction\tn_introns\thas_nearest\t"
        "nearest_gene_id\tnearest_gene_name\tnearest_distance_bp\tnearest_strand\t"
        "nearest_any_gene_id\tnearest_any_gene_name\tnearest_any_distance_bp\tnearest_any_strand\t"
        "locus_mean_depth\tvalley_mean\tvalley_possible\tgap_mean_depth\tn_dup_flag_seen\t"
        "junction_support\tbridge_read_count\tstructure_error\n"
    )
    struct = tmp_path / "s.tsv"
    struct.write_text(
        header + "L1\tc1\tT1\tchr5\t1000\t2000\t+\t2\t400\tGT\tAG\t1.0\t1\ttrue\t"
        "G1\tGm\t40000\t+\tG1\tGm\t40000\t+\t10\t0.0\ttrue\t0.05\t0\t0\t0\n"
        "L1\tt1\tT1\tchr5\t1000\t2000\t+\t2\t400\tGT\tAG\t1.0\t1\ttrue\t"
        "G1\tGm\t40000\t+\tG1\tGm\t40000\t+\t10\t0.0\ttrue\t0.08\t0\t20\t0\n"
        "L1\tt2\tT1\tchr5\t1000\t2000\t+\t2\t400\tGT\tAG\t1.0\t1\ttrue\t"
        "G1\tGm\t40000\t+\tG1\tGm\t40000\t+\t10\t0.0\ttrue\t0.06\t0\t18\t0\n",
        encoding="utf-8",
    )
    tpm = tmp_path / "tpm.tsv"
    tpm.write_text("locus_id\tc1\tc2\tt1\tt2\nL1\t0.0\t0.0\t5.0\t6.0\n")
    counts = tmp_path / "cnt.tsv"
    counts.write_text("locus_id\tc1\tc2\tt1\tt2\nL1\t0\t0\t10\t12\n")
    cfg = TxNovaConfig.model_validate(
        {"genome": {"fasta": "/a.fa", "annotation": "/a.gtf"}, "samples": "/s.tsv"}
    )
    rows = [
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
    cand = apply_gates(
        cfg,
        rows,
        class_tsv=class_tsv,
        reps_tsv=reps,
        structure_tsv=struct,
        locus_tpm=tpm,
        locus_counts=counts,
        out_dir=tmp_path / "cand",
    )
    assert len(cand) == 1
    assert float(cand.iloc[0]["coverage_gap_mean"]) < 0.1


def _four_samples() -> list[SampleRow]:
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


def test_monoexon_passes_without_splice_gate(tmp_path: Path) -> None:
    class_tsv = tmp_path / "c.tsv"
    class_tsv.write_text(
        "transcript_id\tgene_id\tchrom\tstart\tend\tstrand\tn_exons\tlength_nt\texon_structure\tclass\tgene_name\n"
        "T1\tL1\tchr1\t1000\t1400\t+\t1\t400\t1000-1400\tu\tx\n",
        encoding="utf-8",
    )
    reps = tmp_path / "reps.tsv"
    pick_representative_transcript(class_tsv, reps)
    header = (
        "locus_id\tsample_id\ttranscript_id\tchrom\tstart\tend\tstrand\tn_exons\tlength_nt\t"
        "donors\tacceptors\tcanonical_splice_fraction\tn_introns\thas_nearest\t"
        "nearest_gene_id\tnearest_gene_name\tnearest_distance_bp\tnearest_strand\t"
        "nearest_any_gene_id\tnearest_any_gene_name\tnearest_any_distance_bp\tnearest_any_strand\t"
        "locus_mean_depth\tvalley_mean\tvalley_possible\tgap_mean_depth\tn_dup_flag_seen\t"
        "junction_support\tbridge_read_count\tstructure_error\n"
    )
    # n_introns=0 and frac=0.0 is what older structure.tsv wrote.
    struct = tmp_path / "s.tsv"
    struct.write_text(
        header + "L1\tc1\tT1\tchr1\t1000\t1400\t+\t1\t400\t\t\t0.0\t0\tfalse\t"
        "\t\t\t\t\t\t\t\t10\t\t\t\t0\t\t0\n"
        "L1\tt1\tT1\tchr1\t1000\t1400\t+\t1\t400\t\t\t0.0\t0\tfalse\t"
        "\t\t\t\t\t\t\t\t10\t\t\t\t0\t\t0\n"
        "L1\tt2\tT1\tchr1\t1000\t1400\t+\t1\t400\t\t\t0.0\t0\tfalse\t"
        "\t\t\t\t\t\t\t\t10\t\t\t\t0\t\t0\n",
        encoding="utf-8",
    )
    tpm = tmp_path / "tpm.tsv"
    tpm.write_text("locus_id\tc1\tc2\tt1\tt2\nL1\t0.0\t0.0\t5.0\t6.0\n")
    counts = tmp_path / "cnt.tsv"
    counts.write_text("locus_id\tc1\tc2\tt1\tt2\nL1\t0\t0\t10\t12\n")
    cfg = TxNovaConfig.model_validate(
        {"genome": {"fasta": "/a.fa", "annotation": "/a.gtf"}, "samples": "/s.tsv"}
    )
    cand = apply_gates(
        cfg,
        _four_samples(),
        class_tsv=class_tsv,
        reps_tsv=reps,
        structure_tsv=struct,
        locus_tpm=tpm,
        locus_counts=counts,
        out_dir=tmp_path / "cand",
    )
    assert len(cand) == 1
    passed = str(cand.iloc[0]["gates_passed"])
    assert "min_exons" in passed
    assert "canonical_splice" not in passed
    assert int(cand.iloc[0]["n_exons"]) == 1


def test_monoexon_still_fails_distance(tmp_path: Path) -> None:
    class_tsv = tmp_path / "c.tsv"
    class_tsv.write_text(
        "transcript_id\tgene_id\tchrom\tstart\tend\tstrand\tn_exons\tlength_nt\texon_structure\tclass\tgene_name\n"
        "T1\tL1\tchr1\t1000\t1400\t+\t1\t400\t1000-1400\tu\tx\n",
        encoding="utf-8",
    )
    reps = tmp_path / "reps.tsv"
    pick_representative_transcript(class_tsv, reps)
    header = (
        "locus_id\tsample_id\ttranscript_id\tchrom\tstart\tend\tstrand\tn_exons\tlength_nt\t"
        "donors\tacceptors\tcanonical_splice_fraction\tn_introns\thas_nearest\t"
        "nearest_gene_id\tnearest_gene_name\tnearest_distance_bp\tnearest_strand\t"
        "nearest_any_gene_id\tnearest_any_gene_name\tnearest_any_distance_bp\tnearest_any_strand\t"
        "locus_mean_depth\tvalley_mean\tvalley_possible\tgap_mean_depth\tn_dup_flag_seen\t"
        "junction_support\tbridge_read_count\tstructure_error\n"
    )
    struct = tmp_path / "s.tsv"
    struct.write_text(
        header + "L1\tt1\tT1\tchr1\t1000\t1400\t+\t1\t400\t\t\t\t0\ttrue\t"
        "G1\tNear\t800\t+\tG1\tNear\t800\t+\t10\t0.0\ttrue\t0.0\t0\t\t0\n",
        encoding="utf-8",
    )
    tpm = tmp_path / "tpm.tsv"
    tpm.write_text("locus_id\tc1\tc2\tt1\tt2\nL1\t0.0\t0.0\t5.0\t6.0\n")
    counts = tmp_path / "cnt.tsv"
    counts.write_text("locus_id\tc1\tc2\tt1\tt2\nL1\t0\t0\t10\t12\n")
    cfg = TxNovaConfig.model_validate(
        {"genome": {"fasta": "/a.fa", "annotation": "/a.gtf"}, "samples": "/s.tsv"}
    )
    cand = apply_gates(
        cfg,
        _four_samples(),
        class_tsv=class_tsv,
        reps_tsv=reps,
        structure_tsv=struct,
        locus_tpm=tpm,
        locus_counts=counts,
        out_dir=tmp_path / "cand",
    )
    assert cand.empty
