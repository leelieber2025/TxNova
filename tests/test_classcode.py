from __future__ import annotations

from pathlib import Path

import pandas as pd

from txnova.gates import pick_representative_transcript


def _gtf(path: Path, lines: list[str]) -> Path:
    path.write_text("".join(lines), encoding="utf-8")
    return path


def _exon(chrom, start, end, strand, gid, tid, gene_feat=True, gname=None) -> list[str]:
    name = gname or gid
    rows = []
    if gene_feat:
        rows.append(
            f'{chrom}\tX\tgene\t{start}\t{end}\t.\t{strand}\t.\tgene_id "{gid}"; gene_name "{name}";\n'
        )
        rows.append(
            f'{chrom}\tX\ttranscript\t{start}\t{end}\t.\t{strand}\t.\tgene_id "{gid}"; transcript_id "{tid}"; gene_name "{name}";\n'
        )
    rows.append(
        f'{chrom}\tX\texon\t{start}\t{end}\t.\t{strand}\t.\tgene_id "{gid}"; transcript_id "{tid}"; gene_name "{name}";\n'
    )
    return rows


def _classify(tmp: Path, ref_lines: list[str], q_lines: list[str], strand="rf"):
    from txnova import _core

    ref = _gtf(tmp / "ref.gtf", ref_lines)
    mer = _gtf(tmp / "merged.gtf", q_lines)
    out = tmp / "class.tsv"
    r = _core.classify_gtfs(str(mer), str(ref), str(out), f'{{"strandedness":"{strand}"}}')
    df = pd.read_csv(out, sep="\t")
    return r, df


def test_1_pure_intergenic(tmp_path: Path) -> None:
    ref = _exon("chr1", 100, 200, "+", "G1", "T1") + _exon(
        "chr1", 300, 400, "+", "G1", "T1", gene_feat=False
    )
    # fix gene line span — rewrite simpler
    ref = [
        'chr1\tX\tgene\t100\t400\t.\t+\t.\tgene_id "G1"; gene_name "G1";\n',
        'chr1\tX\ttranscript\t100\t400\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n',
        'chr1\tX\texon\t100\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n',
        'chr1\tX\texon\t300\t400\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n',
    ]
    q = [
        'chr1\tX\ttranscript\t1000\t1200\t.\t+\t.\tgene_id "MSTRG.1"; transcript_id "MSTRG.1.1";\n',
        'chr1\tX\texon\t1000\t1080\t.\t+\t.\tgene_id "MSTRG.1"; transcript_id "MSTRG.1.1";\n',
        'chr1\tX\texon\t1120\t1200\t.\t+\t.\tgene_id "MSTRG.1"; transcript_id "MSTRG.1.1";\n',
    ]
    _, df = _classify(tmp_path, ref, q)
    assert df.iloc[0]["class"] == "u"
    assert df.iloc[0]["exon_structure"] == "1000-1080,1120-1200"


def test_2_intron_to_intergenic_is_i(tmp_path: Path) -> None:
    ref = [
        'chr1\tX\tgene\t100\t1000\t.\t+\t.\tgene_id "G1"; gene_name "G1";\n',
        'chr1\tX\ttranscript\t100\t1000\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n',
        'chr1\tX\texon\t100\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n',
        'chr1\tX\texon\t800\t1000\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n',
    ]
    q = [
        'chr1\tX\ttranscript\t400\t1200\t.\t+\t.\tgene_id "MSTRG.1"; transcript_id "MSTRG.1.1";\n',
        'chr1\tX\texon\t400\t450\t.\t+\t.\tgene_id "MSTRG.1"; transcript_id "MSTRG.1.1";\n',
        'chr1\tX\texon\t1100\t1200\t.\t+\t.\tgene_id "MSTRG.1"; transcript_id "MSTRG.1.1";\n',
    ]
    _, df = _classify(tmp_path, ref, q)
    assert df.iloc[0]["class"] == "i"


def test_3_antisense_intronic_is_x(tmp_path: Path) -> None:
    ref = [
        'chr1\tX\tgene\t100\t1000\t.\t+\t.\tgene_id "G1"; gene_name "G1";\n',
        'chr1\tX\ttranscript\t100\t1000\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n',
        'chr1\tX\texon\t100\t150\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n',
        'chr1\tX\texon\t800\t900\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n',
    ]
    q = [
        'chr1\tX\ttranscript\t300\t450\t.\t-\t.\tgene_id "MSTRG.1"; transcript_id "MSTRG.1.1";\n',
        'chr1\tX\texon\t300\t350\t.\t-\t.\tgene_id "MSTRG.1"; transcript_id "MSTRG.1.1";\n',
        'chr1\tX\texon\t400\t450\t.\t-\t.\tgene_id "MSTRG.1"; transcript_id "MSTRG.1.1";\n',
    ]
    _, df = _classify(tmp_path, ref, q)
    assert df.iloc[0]["class"] == "x"


def test_4_dot_strand_overlap_not_u(tmp_path: Path) -> None:
    ref = [
        'chr1\tX\tgene\t100\t200\t.\t+\t.\tgene_id "G1"; gene_name "G1";\n',
        'chr1\tX\ttranscript\t100\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n',
        'chr1\tX\texon\t100\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n',
    ]
    q = [
        'chr1\tX\ttranscript\t150\t250\t.\t.\t.\tgene_id "MSTRG.1"; transcript_id "MSTRG.1.1";\n',
        'chr1\tX\texon\t150\t180\t.\t.\t.\tgene_id "MSTRG.1"; transcript_id "MSTRG.1.1";\n',
        'chr1\tX\texon\t210\t250\t.\t.\t.\tgene_id "MSTRG.1"; transcript_id "MSTRG.1.1";\n',
    ]
    _, df = _classify(tmp_path, ref, q)
    assert df.iloc[0]["class"] != "u"


def test_5_no_gene_line_body_from_exons(tmp_path: Path) -> None:
    ref = [
        'chr1\tX\ttranscript\t100\t400\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n',
        'chr1\tX\texon\t100\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n',
        'chr1\tX\texon\t300\t400\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n',
    ]
    q = [
        'chr1\tX\ttranscript\t250\t280\t.\t+\t.\tgene_id "MSTRG.1"; transcript_id "MSTRG.1.1";\n',
        'chr1\tX\texon\t250\t260\t.\t+\t.\tgene_id "MSTRG.1"; transcript_id "MSTRG.1.1";\n',
        'chr1\tX\texon\t270\t280\t.\t+\t.\tgene_id "MSTRG.1"; transcript_id "MSTRG.1.1";\n',
    ]
    _, df = _classify(tmp_path, ref, q)
    assert df.iloc[0]["class"] == "i"


def test_6_utr_extension(tmp_path: Path) -> None:
    ref = [
        'chr1\tX\tgene\t100\t400\t.\t+\t.\tgene_id "G1"; gene_name "G1";\n',
        'chr1\tX\ttranscript\t100\t400\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n',
        'chr1\tX\texon\t100\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n',
        'chr1\tX\texon\t300\t400\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n',
    ]
    q = [
        'chr1\tX\ttranscript\t300\t500\t.\t+\t.\tgene_id "MSTRG.1"; transcript_id "MSTRG.1.1";\n',
        'chr1\tX\texon\t300\t400\t.\t+\t.\tgene_id "MSTRG.1"; transcript_id "MSTRG.1.1";\n',
        'chr1\tX\texon\t450\t500\t.\t+\t.\tgene_id "MSTRG.1"; transcript_id "MSTRG.1.1";\n',
    ]
    _, df = _classify(tmp_path, ref, q)
    assert df.iloc[0]["class"] == "extension"


def test_7_mixed_locus_excluded(tmp_path: Path) -> None:
    ref = [
        'chr1\tX\tgene\t100\t200\t.\t+\t.\tgene_id "G1"; gene_name "G1";\n',
        'chr1\tX\ttranscript\t100\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n',
        'chr1\tX\texon\t100\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n',
    ]
    q = [
        'chr1\tX\ttranscript\t1000\t1200\t.\t+\t.\tgene_id "MSTRG.1"; transcript_id "MSTRG.1.1";\n',
        'chr1\tX\texon\t1000\t1080\t.\t+\t.\tgene_id "MSTRG.1"; transcript_id "MSTRG.1.1";\n',
        'chr1\tX\texon\t1120\t1200\t.\t+\t.\tgene_id "MSTRG.1"; transcript_id "MSTRG.1.1";\n',
        'chr1\tX\ttranscript\t150\t180\t.\t+\t.\tgene_id "MSTRG.1"; transcript_id "MSTRG.1.2";\n',
        'chr1\tX\texon\t150\t180\t.\t+\t.\tgene_id "MSTRG.1"; transcript_id "MSTRG.1.2";\n',
    ]
    _, df = _classify(tmp_path, ref, q)
    classes = dict(zip(df["transcript_id"], df["class"]))
    assert classes["MSTRG.1.1"] == "u"
    assert classes["MSTRG.1.2"] != "u"
    reps = pick_representative_transcript(tmp_path / "class.tsv", tmp_path / "representatives.tsv")
    assert reps.empty
    loc = pd.read_csv(tmp_path / "locus.class.tsv", sep="\t")
    assert loc.iloc[0]["locus_class"] != "u"
    assert not bool(loc.iloc[0]["all_u"])


def test_bad_cfg_json_fails(tmp_path: Path) -> None:
    from txnova import _core

    ref = tmp_path / "ref.gtf"
    mer = tmp_path / "merged.gtf"
    ref.write_text(
        'chr1\tX\tgene\t100\t200\t.\t+\t.\tgene_id "G1"; gene_name "G1";\n'
        'chr1\tX\ttranscript\t100\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
        'chr1\tX\texon\t100\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n',
        encoding="utf-8",
    )
    mer.write_text(
        'chr1\tX\ttranscript\t1000\t1200\t.\t+\t.\tgene_id "MSTRG.1"; transcript_id "MSTRG.1.1";\n'
        'chr1\tX\texon\t1000\t1080\t.\t+\t.\tgene_id "MSTRG.1"; transcript_id "MSTRG.1.1";\n'
        'chr1\tX\texon\t1120\t1200\t.\t+\t.\tgene_id "MSTRG.1"; transcript_id "MSTRG.1.1";\n',
        encoding="utf-8",
    )
    try:
        _core.classify_gtfs(str(mer), str(ref), str(tmp_path / "c.tsv"), "not-json")
    except Exception as e:
        assert "invalid JSON" in str(e) or "JSON" in str(e)
    else:
        raise AssertionError("expected JSON parse to fail")
