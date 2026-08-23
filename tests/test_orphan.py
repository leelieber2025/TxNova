from __future__ import annotations

from pathlib import Path

import pandas as pd

from txnova.orphan import (
    annotate_orphans,
    collect_orphan_rows,
    is_naming_none,
    parse_exons,
    parse_hmmscan_hits,
    score_conservation,
    summarize_wig,
    ucsc_genome,
)


def test_parse_exons_and_none() -> None:
    assert parse_exons("100-110,200-220") == [(100, 110), (200, 220)]
    assert parse_exons("NA") == []
    assert is_naming_none("none")
    assert is_naming_none("NA")
    assert not is_naming_none("protein_coding")
    assert ucsc_genome("GRCm39") == "mm39"
    assert ucsc_genome("hg38") == "hg38"
    assert ucsc_genome("GRCh38") == "hg38"
    assert ucsc_genome("galGal6") is None


def test_summarize_wig() -> None:
    items = [
        {"start": 10, "end": 11, "value": 0.5},
        {"start": 11, "end": 12, "value": -0.2},
        {"start": 12, "end": 13, "value": 1.0},
    ]
    s = summarize_wig(items)
    assert s["n"] == 3
    assert abs(s["mean"] - 0.433333) < 1e-4
    assert s["max"] == 1.0
    assert abs(s["frac_pos"] - 2 / 3) < 1e-9


def test_parse_hmmscan_hits_keeps_pfam() -> None:
    payload = {
        "status": "SUCCESS",
        "result": {
            "hits": [
                {
                    "name": "MHC_I",
                    "acc": "PF00129",
                    "evalue": 1e-20,
                    "score": 80.1,
                    "desc": "Class I Histocompatibility antigen",
                },
                {
                    "name": "junk",
                    "acc": "PF00000",
                    "evalue": 1.0,
                    "score": 1.0,
                    "desc": "noise",
                },
            ]
        },
    }
    hits = parse_hmmscan_hits(payload)
    assert len(hits) == 1
    assert hits[0]["name"] == "MHC_I"
    assert hits[0]["acc"] == "PF00129"


def test_parse_hmmscan_metadata_name() -> None:
    payload = {
        "status": "SUCCESS",
        "result": {
            "hits": [
                {
                    "name": "000014871",
                    "acc": "PF00129.23",
                    "evalue": 1.4e-18,
                    "score": 67.5,
                    "metadata": {
                        "accession": "PF00129",
                        "identifier": "MHC_I",
                        "description": "Class I Histocompatibility antigen",
                    },
                }
            ]
        },
    }
    hits = parse_hmmscan_hits(payload)
    assert hits[0]["name"] == "MHC_I"
    assert hits[0]["acc"] == "PF00129"
    assert "Histocompatibility" in hits[0]["desc"]


def test_score_conservation_uses_gtf_to_ucsc(tmp_path: Path) -> None:
    seen: list[tuple] = []

    def fake_fetch(genome, track, chrom, start_0, end_0):
        seen.append((genome, track, chrom, start_0, end_0))
        return [{"start": start_0, "end": start_0 + 1, "value": 0.8}]

    out = score_conservation("chr8", [(101, 103)], genome="mm39", fetch=fake_fetch)
    assert seen[0] == ("mm39", "phyloP35way", "chr8", 100, 103)
    assert seen[1][1] == "phastCons35way"
    assert out["phylop_n"] == 1
    assert out["error"] == ""


def test_collect_skips_named() -> None:
    df = pd.DataFrame(
        {
            "locus_id": ["a", "b", "c"],
            "named_overlap": ["none", "protein_coding", "annotated"],
        }
    )
    got = collect_orphan_rows([df])
    assert list(got["locus_id"]) == ["a"]


def test_annotate_orphans_writes_table(tmp_path: Path) -> None:
    src = tmp_path / "unnamed.tsv"
    pd.DataFrame(
        [
            {
                "locus_id": "MSTRG.1",
                "chrom": "chr8",
                "start": 101,
                "end": 200,
                "strand": "+",
                "n_exons": 2,
                "length_nt": 50,
                "exon_structure": "101-110,190-200",
                "named_overlap": "none",
                "longest_orf_aa": 53,
                "coding_label": "hexamer_positive",
            },
            {
                "locus_id": "MSTRG.2",
                "chrom": "chr1",
                "start": 1,
                "end": 10,
                "strand": "+",
                "n_exons": 2,
                "length_nt": 20,
                "exon_structure": "1-5,8-10",
                "named_overlap": "protein_coding",
                "longest_orf_aa": 80,
                "coding_label": "hexamer_positive",
            },
        ]
    ).to_csv(src, sep="\t", index=False)
    fa = tmp_path / "p.fa"
    fa.write_text(">MSTRG.1|x|53\nMELP\n", encoding="utf-8")
    dest = tmp_path / "orphan.tsv"

    def fake_fetch(genome, track, chrom, start_0, end_0):
        return [{"start": start_0, "end": start_0 + 1, "value": 0.2}]

    def fake_scan(seq, *, name="query"):
        assert seq == "MELP"
        assert name == "MSTRG.1"
        return [{"name": "Fake", "acc": "PF0001", "evalue": 1e-10, "score": 40.0, "desc": "toy"}]

    out = annotate_orphans([src], fa, dest, assembly="GRCm39", fetch=fake_fetch, scan=fake_scan)
    assert dest.is_file()
    assert list(out["locus_id"]) == ["MSTRG.1"]
    assert out.iloc[0]["pfam_name"] == "Fake"
    assert out.iloc[0]["pfam_acc"] == "PF0001"
    assert int(out.iloc[0]["phylop_n"]) == 2  # one fake base per exon


def test_annotate_orphans_skips_short_orf(tmp_path: Path) -> None:
    src = tmp_path / "unnamed.tsv"
    pd.DataFrame(
        [
            {
                "locus_id": "MSTRG.1",
                "chrom": "chr8",
                "start": 101,
                "end": 200,
                "strand": "+",
                "n_exons": 1,
                "length_nt": 50,
                "exon_structure": "101-150",
                "named_overlap": "none",
                "longest_orf_aa": 10,
                "coding_label": "hexamer_negative",
            }
        ]
    ).to_csv(src, sep="\t", index=False)
    dest = tmp_path / "orphan.tsv"
    out = annotate_orphans(
        [src],
        None,
        dest,
        assembly="GRCm39",
        min_orf_aa=50,
        fetch=lambda *a, **k: [],
        scan=lambda *a, **k: [],
    )
    assert out.empty
