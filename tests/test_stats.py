from __future__ import annotations

from pathlib import Path

import pandas as pd

from txnova.config import TxNovaConfig
from txnova.samples import SampleRow
from txnova.orchestrator import _join_coding, apply_de_filter, write_final_transcripts
from txnova.stats import de_pass_mask, de_status_series, run_deseq


def _cfg() -> TxNovaConfig:
    return TxNovaConfig.model_validate(
        {"genome": {"fasta": "/a.fa", "annotation": "/a.gtf"}, "samples": "/s.tsv"}
    )


def _rows() -> list[SampleRow]:
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


def test_full_matrix_and_na_not_pass(tmp_path: Path) -> None:
    # High-count known gene + low-count novel (0,0 / 1,1,2,2 is too small — use 0,0 vs 1,2)
    counts = pd.DataFrame(
        {
            "locus_id": ["KNOWN", "NOVEL", "ZERO"],
            "c1": [1000, 0, 0],
            "c2": [1100, 0, 0],
            "t1": [1000, 1, 0],
            "t2": [900, 2, 0],
        }
    )
    path = tmp_path / "locus_counts.tsv"
    counts.to_csv(path, sep="\t", index=False)
    res = run_deseq(_cfg(), _rows(), path, tmp_path / "de.tsv")
    assert set(res["locus_id"]) == {"KNOWN", "NOVEL", "ZERO"}
    zero = res.set_index("locus_id").loc["ZERO"]
    assert pd.isna(zero["padj"])
    mask = de_pass_mask(res, _cfg())
    status = de_status_series(res, _cfg())
    # all-zero: not fitted — padj and pvalue NA — must not pass
    assert not bool(mask[res["locus_id"] == "ZERO"].iloc[0])
    assert status[res["locus_id"] == "ZERO"].iloc[0] == ""
    novel = res.set_index("locus_id").loc["NOVEL"]
    novel_pass = bool(mask[res["locus_id"] == "NOVEL"].iloc[0])
    novel_status = status[res["locus_id"] == "NOVEL"].iloc[0]
    if (
        pd.notna(novel["padj"])
        and float(novel["padj"]) < 0.05
        and float(novel["log2FoldChange"]) >= 2.0
    ):
        assert novel_pass
        assert novel_status == "wald"
    elif (
        pd.isna(novel["padj"])
        and pd.notna(novel["pvalue"])
        and float(novel["log2FoldChange"]) >= 2.0
    ):
        assert novel_pass
        assert novel_status == "low_count"
    else:
        assert not novel_pass
        assert novel_status == ""


def test_de_filter_from_gates_view_is_repeatable(tmp_path: Path) -> None:
    gates = tmp_path / "candidates.gates.tsv"
    gates.write_text(
        "locus_id\tlocus_coord\tchrom\tstart\tend\tstrand\tn_exons\tlength_nt\t"
        "exon_structure\tclass\tgates_passed\trepresentative_transcript_id\n"
        "L1\tchr1:1-10:+\tchr1\t1\t10\t+\t2\t400\t1-10\tu\tclass\tT1\n",
        encoding="utf-8",
    )
    res = pd.DataFrame(
        {"locus_id": ["L1"], "log2FoldChange": [3.0], "pvalue": [1e-5], "padj": [0.001]}
    )
    dest = tmp_path / "candidates.de.tsv"
    apply_de_filter(_cfg(), _rows(), gates, res, dest)
    apply_de_filter(_cfg(), _rows(), gates, res, dest)
    out = pd.read_csv(dest, sep="\t")
    assert "padj" in out.columns
    assert "padj_x" not in out.columns
    assert list(out["locus_id"]) == ["L1"]
    assert list(out["de_status"]) == ["wald"]


def test_low_count_independent_filter_passes(tmp_path: Path) -> None:
    gates = tmp_path / "candidates.gates.tsv"
    gates.write_text(
        "locus_id\tlocus_coord\tchrom\tstart\tend\tstrand\tn_exons\tlength_nt\t"
        "exon_structure\tclass\tgates_passed\trepresentative_transcript_id\n"
        "LOW\tchr1:1-10:+\tchr1\t1\t10\t+\t2\t400\t1-10\tu\tclass,detect\tT1\n"
        "COOK\tchr1:20-30:+\tchr1\t20\t30\t+\t2\t400\t20-30\tu\tclass,detect\tT2\n"
        "FAIL\tchr1:40-50:+\tchr1\t40\t50\t+\t2\t400\t40-50\tu\tclass,detect\tT3\n",
        encoding="utf-8",
    )
    res = pd.DataFrame(
        {
            "locus_id": ["LOW", "COOK", "FAIL"],
            "baseMean": [0.8, 50.0, 12.0],
            "log2FoldChange": [3.5, 4.0, 2.5],
            "pvalue": [0.14, float("nan"), 0.20],
            "padj": [float("nan"), float("nan"), 0.40],
        }
    )
    dest = tmp_path / "candidates.de.tsv"
    apply_de_filter(_cfg(), _rows(), gates, res, dest)
    out = pd.read_csv(dest, sep="\t")
    assert list(out["locus_id"]) == ["LOW"]
    assert list(out["de_status"]) == ["low_count"]
    assert bool(out.iloc[0]["de_pass"])


def test_de_status_constructed_rows() -> None:
    cfg = _cfg()
    res = pd.DataFrame(
        {
            "locus_id": ["W", "L", "Z", "NS", "DOWN"],
            "log2FoldChange": [3.0, 0.6, float("nan"), 3.0, 0.4],
            "pvalue": [1e-6, 0.14, float("nan"), 0.02, 0.3],
            "padj": [1e-4, float("nan"), float("nan"), 0.20, float("nan")],
        }
    )
    status = de_status_series(res, cfg)
    mask = de_pass_mask(res, cfg)
    assert list(status) == ["wald", "low_count", "", "", ""]
    assert list(mask) == [True, True, False, False, False]


def test_de_status_requires_pvalue() -> None:
    from txnova.errors import TxNovaError

    res = pd.DataFrame({"locus_id": ["W"], "log2FoldChange": [3.0], "padj": [1e-4]})
    try:
        de_status_series(res, _cfg())
    except TxNovaError as e:
        assert "pvalue" in str(e)
    else:
        raise AssertionError("expected TxNovaError")


def test_join_coding_twice_no_suffix_columns(tmp_path: Path) -> None:
    src = tmp_path / "candidates.de.tsv"
    src.write_text(
        "locus_id\tlocus_coord\tchrom\tstart\tend\tstrand\tn_exons\tlength_nt\t"
        "exon_structure\tclass\tgates_passed\trepresentative_transcript_id\n"
        "L1\tchr1:1-10:+\tchr1\t1\t10\t+\t2\t400\t1-10\tu\tclass\tT1\n",
        encoding="utf-8",
    )
    orfs = tmp_path / "orfs.tsv"
    orfs.write_text(
        "locus_id\ttranscript_id\tlongest_orf_aa\torf_complete\tcoding_score\tfickett_score\tcoding_label\n"
        "L1\tT1\t80\ttrue\t0.2\t0.5\tcoding\n",
        encoding="utf-8",
    )
    dest = tmp_path / "candidates.tsv"
    cfg = _cfg()
    _join_coding(cfg, src, orfs, dest, _rows())
    first = dest.read_text()
    _join_coding(cfg, src, orfs, dest, _rows())
    assert dest.read_text() == first
    _join_coding(cfg, dest, orfs, dest, _rows())
    out = pd.read_csv(dest, sep="\t")
    assert "coding_label" in out.columns
    assert "coding_label_x" not in out.columns
    assert str(out.iloc[0]["coding_label"]) == "coding"


def test_final_transcripts_follow_candidates_not_gates(tmp_path: Path) -> None:
    cls = tmp_path / "class.tsv"
    cls.write_text(
        "transcript_id\tgene_id\tclass\nT1\tL1\tu\nT2\tL2\tu\n",
        encoding="utf-8",
    )
    cand = tmp_path / "candidates.tsv"
    cand.write_text("locus_id\nL1\n", encoding="utf-8")
    dest = tmp_path / "transcripts.tsv"
    write_final_transcripts(cls, cand, dest)
    out = pd.read_csv(dest, sep="\t")
    assert list(out["gene_id"]) == ["L1"]
