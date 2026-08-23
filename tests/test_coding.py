from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from txnova.config import packaged_hexamer_table


def _write_fasta(path: Path, chrom: str, seq: str) -> None:
    width = 60
    chunks = [seq[i : i + width] for i in range(0, len(seq), width)]
    header = f">{chrom}\n"
    body = "".join(c + "\n" for c in chunks)
    path.write_text(header + body, encoding="utf-8")
    offset = len(header.encode())
    linebases = width
    linewidth = width + 1
    path.with_suffix(path.suffix + ".fai").write_text(
        f"{chrom}\t{len(seq)}\t{offset}\t{linebases}\t{linewidth}\n",
        encoding="utf-8",
    )


def _gtf(chrom: str, tid: str, gid: str, end: int) -> str:
    return (
        f'{chrom}\tX\ttranscript\t1\t{end}\t.\t+\t.\tgene_id "{gid}"; transcript_id "{tid}";\n'
        f'{chrom}\tX\texon\t1\t{end}\t.\t+\t.\tgene_id "{gid}"; transcript_id "{tid}";\n'
    )


def test_hexamer_llr_known_vs_noncoding(tmp_path: Path) -> None:
    from txnova import _core

    table = packaged_hexamer_table()
    coding_seq = "ATG" + ("GAACGG" * 30) + "TAA"
    nc_seq = "ATG" + ("CTTCTT" * 30) + "TAA"
    pad = "N" * 20
    chrom_seq = coding_seq + pad + nc_seq
    fa = tmp_path / "g.fa"
    _write_fasta(fa, "chr1", chrom_seq)
    start_nc = len(coding_seq) + len(pad) + 1
    end_nc = start_nc + len(nc_seq) - 1
    gtf = tmp_path / "m.gtf"
    gtf.write_text(
        _gtf("chr1", "Tcod", "Lcod", len(coding_seq))
        + (
            f"chr1\tX\ttranscript\t{start_nc}\t{end_nc}\t.\t+\t.\t"
            f'gene_id "Lnc"; transcript_id "Tnc";\n'
            f"chr1\tX\texon\t{start_nc}\t{end_nc}\t.\t+\t.\t"
            f'gene_id "Lnc"; transcript_id "Tnc";\n'
        ),
        encoding="utf-8",
    )
    reps = tmp_path / "reps.tsv"
    reps.write_text("locus_id\ttranscript_id\nLcod\tTcod\nLnc\tTnc\n", encoding="utf-8")
    orfs = tmp_path / "orfs.tsv"
    peps = tmp_path / "pep.fa"
    n = _core.scan_orfs(
        str(fa),
        str(gtf),
        str(reps),
        str(orfs),
        str(peps),
        json.dumps(
            {
                "min_orf_aa": 50,
                "hexamer_table": str(table),
                "hexamer_coding_min": 0.0,
                "hexamer_noncoding_max": 0.0,
            }
        ),
    )
    assert n["n_orfs"] == 2
    df = pd.read_csv(orfs, sep="\t")
    by = df.set_index("locus_id")
    assert int(by.loc["Lcod", "longest_orf_aa"]) >= 50
    assert str(by.loc["Lcod", "orf_complete"]).lower() == "true"
    assert float(by.loc["Lcod", "coding_score"]) > 0
    assert by.loc["Lcod", "coding_label"] == "hexamer_positive"
    assert float(by.loc["Lnc", "coding_score"]) < 0
    assert by.loc["Lnc", "coding_label"] == "hexamer_negative"
    assert pd.notna(by.loc["Lcod", "fickett_score"])


def test_no_orf_is_noncoding(tmp_path: Path) -> None:
    from txnova import _core

    seq = "CCCC" * 80
    fa = tmp_path / "g.fa"
    _write_fasta(fa, "chr1", seq)
    gtf = tmp_path / "m.gtf"
    gtf.write_text(_gtf("chr1", "T1", "L1", len(seq)), encoding="utf-8")
    reps = tmp_path / "reps.tsv"
    reps.write_text("locus_id\ttranscript_id\nL1\tT1\n", encoding="utf-8")
    orfs = tmp_path / "orfs.tsv"
    _core.scan_orfs(
        str(fa),
        str(gtf),
        str(reps),
        str(orfs),
        str(tmp_path / "pep.fa"),
        json.dumps(
            {
                "min_orf_aa": 50,
                "hexamer_table": str(packaged_hexamer_table()),
                "hexamer_coding_min": 0.0,
                "hexamer_noncoding_max": 0.0,
            }
        ),
    )
    df = pd.read_csv(orfs, sep="\t")
    assert int(df.iloc[0]["longest_orf_aa"]) == 0
    assert df.iloc[0]["coding_score"] == "NA" or pd.isna(df.iloc[0]["coding_score"])
    assert df.iloc[0]["coding_label"] == "no_orf"


def test_short_orf_is_reported_not_zero(tmp_path: Path) -> None:
    from txnova import _core

    # 20 aa complete ORF — below the old min_orf_aa=50 reporting floor.
    seq = "ATG" + ("GAA" * 19) + "TAA" + "C" * 40
    fa = tmp_path / "g.fa"
    _write_fasta(fa, "chr1", seq)
    gtf = tmp_path / "m.gtf"
    gtf.write_text(_gtf("chr1", "T1", "L1", len(seq)), encoding="utf-8")
    reps = tmp_path / "reps.tsv"
    reps.write_text("locus_id\ttranscript_id\nL1\tT1\n", encoding="utf-8")
    orfs = tmp_path / "orfs.tsv"
    _core.scan_orfs(
        str(fa),
        str(gtf),
        str(reps),
        str(orfs),
        str(tmp_path / "pep.fa"),
        json.dumps(
            {
                "min_orf_aa": 50,
                "hexamer_table": str(packaged_hexamer_table()),
                "hexamer_coding_min": 0.0,
                "hexamer_noncoding_max": 0.0,
            }
        ),
    )
    df = pd.read_csv(orfs, sep="\t")
    assert int(df.iloc[0]["longest_orf_aa"]) == 20
    assert pd.notna(df.iloc[0]["coding_score"])
