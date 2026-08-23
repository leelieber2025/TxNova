from __future__ import annotations

from pathlib import Path

import pandas as pd

from txnova.gene_score import (
    chrom_class,
    junction_min,
    rank_loci,
    stage2_adjust,
    top_rank_ids,
    write_gene_rank,
)


def test_chrom_and_junction_helpers() -> None:
    assert chrom_class("chr11") == "primary"
    assert chrom_class("GL456233.2") == "unplaced"
    assert chrom_class("MU069435.1") == "unplaced"
    assert chrom_class("chrY") == "primary"
    assert chrom_class("chr22") == "primary"
    assert junction_min("25,22,0,16") == 0
    assert junction_min("783,825") == 783
    assert junction_min("") is None


def test_top_rank_ids(tmp_path: Path) -> None:
    p = tmp_path / "gene_rank.tsv"
    p.write_text(
        "gene_rank\tlocus_id\n2\tB\n1\tA\n3\tC\n",
        encoding="utf-8",
    )
    assert top_rank_ids(p, 2) == ["A", "B"]


def test_unplaced_and_nosplice_rank_below_chromosomal() -> None:
    rows = [
        {
            "locus_id": "U",
            "chrom": "GL456233.2",
            "start": 1,
            "end": 20000,
            "n_exons": 8,
            "length_nt": 1300,
            "junction_support": "20,20,20,20,20,20,20",
            "longest_orf_aa": 200,
            "orf_complete": True,
            "coding_score": 0.1,
            "coding_label": "hexamer_positive",
            "nearest_distance_bp": 1e6,
            "treat_n_detected": 6,
        },
        {
            "locus_id": "C",
            "chrom": "chr3",
            "start": 100,
            "end": 1800,
            "n_exons": 1,
            "length_nt": 1682,
            "junction_support": "",
            "longest_orf_aa": 205,
            "orf_complete": True,
            "coding_score": 0.33,
            "coding_label": "hexamer_positive",
            "nearest_distance_bp": 3900,
            "nearest_any_distance_bp": 3900,
            "treat_n_detected": 6,
        },
        {
            "locus_id": "Z",
            "chrom": "chr4",
            "start": 1,
            "end": 2000,
            "n_exons": 2,
            "length_nt": 500,
            "junction_support": "0",
            "longest_orf_aa": 0,
            "coding_score": 0,
            "coding_label": "hexamer_negative",
            "nearest_distance_bp": 20000,
            "treat_n_detected": 6,
        },
    ]
    out = rank_loci(pd.DataFrame(rows))
    assert list(out["locus_id"])[0] == "C"
    assert "unplaced" in str(out.set_index("locus_id").loc["U", "gene_flags"])
    assert "no_splice" in str(out.set_index("locus_id").loc["Z", "gene_flags"])


def test_stage2_te_and_processed_penalties() -> None:
    rec = {"n_exons": 1, "longest_orf_aa": 80}
    te, fl = stage2_adjust(rec, pfam_name="Gag_p30", pfam_desc="Gag P30", phylop_mean=0.1)
    assert te <= -50
    assert "te" in fl
    pr, fl2 = stage2_adjust(
        rec, pfam_name="Mt_ATP-synt_D", pfam_desc="ATP synthase D", phylop_mean=0.9
    )
    assert pr <= -30
    assert "processed" in fl2


def test_write_gene_rank_uses_function_and_resorts(tmp_path: Path) -> None:
    src = tmp_path / "u.tsv"
    pd.DataFrame(
        [
            {
                "locus_id": "TE",
                "locus_coord": "chr1:1-400:+",
                "chrom": "chr1",
                "start": 1,
                "end": 400,
                "n_exons": 1,
                "length_nt": 400,
                "exon_structure": "1-400",
                "nearest_gene_name": "X",
                "nearest_distance_bp": 20000,
                "treat_median_tpm": 10,
                "control_max_tpm": 10,
                "treat_n_detected": 6,
                "junction_support": "",
                "longest_orf_aa": 90,
                "orf_complete": True,
                "coding_score": 0.5,
                "coding_label": "hexamer_positive",
            },
            {
                "locus_id": "OK",
                "locus_coord": "chr2:1-800:+",
                "chrom": "chr2",
                "start": 1,
                "end": 800,
                "n_exons": 1,
                "length_nt": 800,
                "exon_structure": "1-800",
                "nearest_gene_name": "Y",
                "nearest_distance_bp": 20000,
                "treat_median_tpm": 2,
                "control_max_tpm": 2,
                "treat_n_detected": 6,
                "junction_support": "",
                "longest_orf_aa": 120,
                "orf_complete": True,
                "coding_score": 0.3,
                "coding_label": "hexamer_positive",
            },
        ]
    ).to_csv(src, sep="\t", index=False)

    def fake_scan(seq: str, name: str = "q") -> list[dict]:
        if name == "TE":
            return [
                {"name": "Gag_p30", "acc": "PF02093", "evalue": 1e-20, "score": 80, "desc": "Gag"}
            ]
        return []

    def fake_fetch(*_a, **_k) -> list[dict]:
        return [{"start": 0, "end": 1, "value": 0.1}]

    peps = tmp_path / "pep.fa"
    peps.write_text(">TE|t|90\nM" + "A" * 89 + "\n>OK|o|120\nM" + "G" * 119 + "\n")
    dest = tmp_path / "rank.tsv"
    out = write_gene_rank(
        [src],
        dest,
        peptides_fa=peps,
        assembly="GRCm39",
        scan=fake_scan,
        fetch=fake_fetch,
    )
    assert dest.is_file()
    assert list(out["locus_id"])[0] == "OK"
    te = out.set_index("locus_id").loc["TE"]
    assert "te" in str(te["gene_flags"])
    assert te["pfam_name"] == "Gag_p30"


def test_stage2_can_write_float_onto_initially_empty_column() -> None:
    ranked = pd.DataFrame({"locus_id": ["A"], "phylop_mean": [""]})
    ranked["phylop_mean"] = ranked["phylop_mean"].astype(object)
    ranked.loc[ranked["locus_id"] == "A", "phylop_mean"] = 1.87
    assert float(ranked.iloc[0]["phylop_mean"]) == 1.87
