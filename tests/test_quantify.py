from __future__ import annotations

from pathlib import Path

import pandas as pd

from txnova import _core

FIXTURES = Path(__file__).resolve().parent / "fixtures"


def _gtf(path: Path, text: str) -> Path:
    path.write_text(text)
    return path


def _quant(tmp: Path, mer: Path, bam: Path) -> Path:
    samples = '{"samples":[{"sample_id":"s1","bam":"%s","dup_flag_seen":false}]}' % bam
    out = tmp / "q"
    _core.quantify_gtf(str(mer), samples, str(out), '{"strandedness":"unstranded","min_mapq":0}')
    return out


def test_intron_only_not_counted(tmp_path: Path) -> None:
    mer = _gtf(
        tmp_path / "m.gtf",
        'chr1\tX\ttranscript\t100\t400\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
        'chr1\tX\texon\t100\t150\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
        'chr1\tX\texon\t350\t400\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n',
    )
    out = _quant(tmp_path, mer, FIXTURES / "quant_intron_only.bam")
    loc = pd.read_csv(out / "locus_counts.tsv", sep="\t")
    assert float(loc.iloc[0]["s1"]) == 0.0


def test_exon_hit_counted(tmp_path: Path) -> None:
    mer = _gtf(
        tmp_path / "m.gtf",
        'chr1\tX\ttranscript\t100\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
        'chr1\tX\texon\t100\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n',
    )
    out = _quant(tmp_path, mer, FIXTURES / "quant_exon_hit.bam")
    loc = pd.read_csv(out / "locus_counts.tsv", sep="\t")
    assert float(loc.iloc[0]["s1"]) == 1.0


def test_cross_locus_discard(tmp_path: Path) -> None:
    mer = _gtf(
        tmp_path / "m.gtf",
        'chr1\tX\ttranscript\t100\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
        'chr1\tX\texon\t100\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
        'chr1\tX\ttranscript\t150\t250\t.\t+\t.\tgene_id "G2"; transcript_id "T2";\n'
        'chr1\tX\texon\t150\t250\t.\t+\t.\tgene_id "G2"; transcript_id "T2";\n',
    )
    out = _quant(tmp_path, mer, FIXTURES / "quant_cross_locus.bam")
    loc = pd.read_csv(out / "locus_counts.tsv", sep="\t")
    assert loc["s1"].astype(float).sum() == 0.0


def test_tpm_denominator_full_table(tmp_path: Path) -> None:
    mer = _gtf(
        tmp_path / "m.gtf",
        'chr1\tX\ttranscript\t100\t200\t.\t+\t.\tgene_id "KNOWN"; transcript_id "TK";\n'
        'chr1\tX\texon\t100\t200\t.\t+\t.\tgene_id "KNOWN"; transcript_id "TK";\n'
        'chr1\tX\ttranscript\t500\t600\t.\t+\t.\tgene_id "NOVEL"; transcript_id "TN";\n'
        'chr1\tX\texon\t500\t600\t.\t+\t.\tgene_id "NOVEL"; transcript_id "TN";\n',
    )
    out = _quant(tmp_path, mer, FIXTURES / "quant_tpm_novel.bam")
    tpm1 = pd.read_csv(out / "locus_tpm.tsv", sep="\t").set_index("locus_id")
    novel_only = float(tpm1.loc["NOVEL", "s1"])
    samples2 = '{"samples":[{"sample_id":"s1","bam":"%s","dup_flag_seen":false}]}' % (
        FIXTURES / "quant_tpm_both.bam"
    )
    out2 = tmp_path / "q2"
    _core.quantify_gtf(str(mer), samples2, str(out2), '{"strandedness":"unstranded","min_mapq":0}')
    tpm2 = pd.read_csv(out2 / "locus_tpm.tsv", sep="\t").set_index("locus_id")
    assert float(tpm2.loc["NOVEL", "s1"]) < novel_only
    assert "KNOWN" in tpm2.index


def test_secondary_ignored(tmp_path: Path) -> None:
    mer = _gtf(
        tmp_path / "m.gtf",
        'chr1\tX\ttranscript\t100\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
        'chr1\tX\texon\t100\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n',
    )
    out = _quant(tmp_path, mer, FIXTURES / "quant_secondary.bam")
    loc = pd.read_csv(out / "locus_counts.tsv", sep="\t")
    assert float(loc.iloc[0]["s1"]) == 0.0


def test_bad_cfg_json_fails(tmp_path: Path) -> None:
    mer = _gtf(
        tmp_path / "m.gtf",
        'chr1\tX\ttranscript\t100\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
        'chr1\tX\texon\t100\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n',
    )
    samples = '{"samples":[{"sample_id":"s1","bam":"%s","dup_flag_seen":false}]}' % (
        FIXTURES / "quant_exon_hit.bam"
    )
    try:
        _core.quantify_gtf(str(mer), samples, str(tmp_path / "q"), "not-json")
    except Exception as e:
        assert "invalid JSON" in str(e) or "JSON" in str(e)
    else:
        raise AssertionError("expected JSON parse to fail")
