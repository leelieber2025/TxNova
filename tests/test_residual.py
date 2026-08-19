from __future__ import annotations

import pandas as pd

from pathlib import Path

from txnova.errors import TxNovaError
from txnova.residual import (
    apply_terminal_extents,
    clip_terminals_from_genes,
    clip_to_contig,
    close_to_small_rna,
    cluster_leak,
    cluster_junctions,
    exon_structure,
    merge_gene_bodies,
    overlaps_gene_body,
    same_locus,
    select_junctions,
    write_residual_gtf,
    write_residuals,
)


def test_shared_site_and_chain() -> None:
    a = {"chrom": "chr1", "start": 100, "end": 200, "strand": "+"}
    b = {"chrom": "chr1", "start": 100, "end": 300, "strand": "+"}  # shared donor
    c = {"chrom": "chr1", "start": 250, "end": 400, "strand": "+"}  # chain after a: gap=49
    d = {"chrom": "chr1", "start": 50000, "end": 50100, "strand": "+"}  # beyond MAX_EXON
    assert same_locus(a, b)
    assert same_locus(a, c)
    assert not same_locus(a, d)


def test_cluster_chains_three_exons() -> None:
    rows = [
        {
            "chrom": "chr8",
            "start": 1000,
            "end": 1100,
            "strand": "-",
            "treat_sum": 8,
            "treat_n_detected": 4,
            "control_max": 0,
            "nearest_distance_bp": 20000,
            "nearest_gene_name": "X",
        },
        {
            "chrom": "chr8",
            "start": 1200,
            "end": 1300,
            "strand": "-",
            "treat_sum": 6,
            "treat_n_detected": 3,
            "control_max": 0,
            "nearest_distance_bp": 19800,
            "nearest_gene_name": "X",
        },
    ]
    groups = cluster_junctions(rows)
    assert len(groups) == 1
    assert len(groups[0]) == 2
    struct, n_exons, s, e = exon_structure(rows)
    assert n_exons == 3
    assert "1101-1199" in struct


def test_cluster_leak_drops_in_gene_keeps_800bp() -> None:
    leak = pd.DataFrame(
        [
            {
                "chrom": "chr1",
                "start": 100,
                "end": 200,
                "strand": "+",
                "status": "unassembled",
                "overlaps_gene": False,
                "treat_sum": 5,
                "treat_n_detected": 3,
                "control_max": 0,
                "nearest_distance_bp": 800,
                "nearest_gene_name": "Near",
            },
            {
                "chrom": "chr2",
                "start": 500,
                "end": 600,
                "strand": "+",
                "status": "unassembled",
                "overlaps_gene": True,
                "treat_sum": 20,
                "treat_n_detected": 6,
                "control_max": 0,
                "nearest_distance_bp": 0,
                "nearest_gene_name": "InGene",
            },
            {
                "chrom": "chr3",
                "start": 1000,
                "end": 1100,
                "strand": "+",
                "status": "unassembled",
                "overlaps_gene": False,
                "treat_sum": 9,
                "treat_n_detected": 4,
                "control_max": 0,
                "nearest_distance_bp": 12000,
                "nearest_gene_name": "Far",
            },
        ]
    )
    out, _ = cluster_leak(leak)
    # overlaps_gene dropped. 800 bp is a gate, not a harvest cut.
    assert list(out["chrom"]) == ["chr3", "chr1"]
    assert int(out.iloc[0]["n_exons"]) == 2
    assert out.iloc[0]["residual_id"] == "RSDL.1"


def test_write_residuals_empty_and_populated(tmp_path: Path) -> None:
    dest = tmp_path / "residual.tsv"
    empty = write_residuals(tmp_path / "missing.tsv", dest)
    assert empty.empty
    assert dest.is_file()
    header = dest.read_text(encoding="utf-8").splitlines()[0].split("\t")
    assert "residual_id" in header
    assert "n_shared_site" in header
    src = tmp_path / "leak.tsv"
    pd.DataFrame(
        [
            {
                "chrom": "chr4",
                "start": 1000,
                "end": 1100,
                "strand": "+",
                "status": "unassembled",
                "overlaps_gene": False,
                "treat_sum": 8,
                "treat_n_detected": 3,
                "control_max": 0,
                "nearest_distance_bp": 12000,
                "nearest_gene_name": "Far",
            },
            {
                "chrom": "chr4",
                "start": 1200,
                "end": 1300,
                "strand": "+",
                "status": "unassembled",
                "overlaps_gene": False,
                "treat_sum": 6,
                "treat_n_detected": 3,
                "control_max": 0,
                "nearest_distance_bp": 11800,
                "nearest_gene_name": "Far",
            },
        ]
    ).to_csv(src, sep="\t", index=False)
    out = write_residuals(src, dest)
    assert len(out) == 1
    assert int(out.iloc[0]["n_junctions"]) == 2
    assert int(out.iloc[0]["n_exons"]) == 3
    assert "1101-1199" in str(out.iloc[0]["exon_structure"])
    written = pd.read_csv(dest, sep="\t")
    assert list(written["residual_id"]) == ["RSDL.1"]
    gtf = dest.with_name("residual.gtf")
    assert gtf.is_file()
    text = gtf.read_text()
    assert 'gene_id "RSDL.1"' in text
    assert "txnova\texon\t" in text


def test_universe_concat(tmp_path: Path) -> None:
    from txnova.residual import write_universe_gtf

    merged = tmp_path / "merged.gtf"
    merged.write_text('chr1\tX\tgene\t1\t10\t.\t+\t.\tgene_id "MSTRG.1";\n', encoding="utf-8")
    res = tmp_path / "residual.gtf"
    res.write_text('chr2\ttxnova\tgene\t5\t8\t.\t-\t.\tgene_id "RSDL.1";\n', encoding="utf-8")
    dest = tmp_path / "universe.gtf"
    write_universe_gtf(merged, res, dest)
    body = dest.read_text()
    assert "MSTRG.1" in body
    assert "RSDL.1" in body


def test_universe_drops_cds_utr(tmp_path: Path) -> None:
    from txnova.residual import write_universe_gtf

    merged = tmp_path / "merged.gtf"
    merged.write_text(
        'chr1\tX\tgene\t1\t20\t.\t+\t.\tgene_id "G";\n'
        'chr1\tX\ttranscript\t1\t20\t.\t+\t.\tgene_id "G"; transcript_id "T";\n'
        'chr1\tX\texon\t1\t20\t.\t+\t.\tgene_id "G"; transcript_id "T";\n'
        'chr1\tX\tCDS\t5\t16\t.\t+\t.\tgene_id "G"; transcript_id "T";\n'
        'chr1\tX\tUTR\t1\t4\t.\t+\t.\tgene_id "G"; transcript_id "T";\n',
        encoding="utf-8",
    )
    res = tmp_path / "residual.gtf"
    res.write_text("", encoding="utf-8")
    dest = tmp_path / "universe.gtf"
    write_universe_gtf(merged, res, dest)
    body = dest.read_text()
    assert "\tCDS\t" not in body
    assert "\tUTR\t" not in body
    assert "\texon\t" in body


def _junc(chrom, start, end, strand="+", **kw):
    row = {
        "chrom": chrom,
        "start": start,
        "end": end,
        "strand": strand,
        "status": "unassembled",
        "overlaps_gene": False,
        "treat_sum": 8,
        "treat_n_detected": 3,
        "control_max": 0,
        "nearest_distance_bp": 12000,
        "nearest_gene_name": "Far",
    }
    row.update(kw)
    return row


def test_drops_antisense_gene_body_keeps_nearby_nonoverlap(tmp_path: Path) -> None:
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(
        'chr5\tX\tgene\t1000\t2000\t.\t+\t.\tgene_id "G1"; gene_name "Col"; gene_type "protein_coding";\n'
        'chr6\tX\tgene\t5000\t6000\t.\t-\t.\tgene_id "G2"; gene_name "NearOpp"; gene_type "protein_coding";\n',
        encoding="utf-8",
    )
    genes = merge_gene_bodies([gtf])
    leak = pd.DataFrame(
        [
            _junc("chr5", 1200, 1800, "-"),  # inside Col, opposite strand
            _junc("chr6", 3000, 3100, "+"),  # 1899 bp from opposite-strand gene, no overlap
        ]
    )
    out, _ = cluster_leak(leak, genes=genes)
    assert list(out["chrom"]) == ["chr6"]


def test_does_not_chain_across_a_remaining_gene(tmp_path: Path) -> None:
    """A gene in the gap blocks chaining; each side stays its own locus.

    Previously the two junctions merged, the span overlapped the gene, and
    both sides were discarded.
    """
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(
        'chr7\tX\tgene\t2000\t4000\t.\t+\t.\tgene_id "G"; gene_name "Block"; gene_type "protein_coding";\n',
        encoding="utf-8",
    )
    leak = pd.DataFrame(
        [
            _junc("chr7", 1000, 1100, "-"),
            _junc("chr7", 5000, 5100, "-"),
        ]
    )
    out, _ = cluster_leak(leak, genes=merge_gene_bodies([gtf]))
    assert len(out) == 2
    starts = sorted(int(x) for x in out["start"])
    # 30 nt stubs: 1000-30=970 and 5000-30=4970
    assert starts == [970, 4970]


def test_naming_gtf_does_not_apply_distance_either_strand(tmp_path: Path) -> None:
    """Close opposite-strand gene is not a harvest distance gate."""
    gtf = tmp_path / "m39.gtf"
    gtf.write_text(
        'chr8\tX\tgene\t4000\t4500\t.\t-\t.\tgene_id "Gm"; gene_name "GmClose"; gene_type "lncRNA";\n',
        encoding="utf-8",
    )
    leak = pd.DataFrame([_junc("chr8", 1000, 1100, "+", nearest_distance_bp=20000)])
    out, _ = cluster_leak(leak, genes=merge_gene_bodies([gtf]))
    assert len(out) == 1
    assert out.iloc[0]["chrom"] == "chr8"


def test_drops_same_strand_stuck_to_gene_keeps_1kb_and_opposite(tmp_path: Path) -> None:
    gtf = tmp_path / "m39.gtf"
    gtf.write_text(
        'chr9\tX\tgene\t2000\t3000\t.\t+\t.\tgene_id "A"; gene_name "Akap"; gene_type "protein_coding";\n'
        'chr10\tX\tgene\t5000\t6000\t.\t+\t.\tgene_id "B"; gene_name "FarSame"; gene_type "lncRNA";\n'
        'chr11\tX\tgene\t1120\t2000\t.\t-\t.\tgene_id "C"; gene_name "Opp17"; gene_type "protein_coding";\n',
        encoding="utf-8",
    )
    genes = merge_gene_bodies([gtf])
    leak = pd.DataFrame(
        [
            _junc("chr9", 1000, 1980, "+"),  # 19 bp upstream of Akap, same strand
            _junc("chr10", 1000, 1100, "+"),  # 3899 bp from FarSame
            _junc("chr11", 1000, 1100, "+"),  # 19 bp from opposite-strand gene
        ]
    )
    out, _ = cluster_leak(leak, genes=genes)
    assert set(out["chrom"]) == {"chr10", "chr11"}


def test_harvest_keeps_800bp_same_strand(tmp_path: Path) -> None:
    """200 nt knife only. 5 kb is a gate, not a harvest cut."""
    gtf = tmp_path / "m39.gtf"
    gtf.write_text(
        'chr1\tX\tgene\t2000\t3000\t.\t+\t.\tgene_id "G"; gene_name "Near"; gene_type "protein_coding";\n',
        encoding="utf-8",
    )
    leak = pd.DataFrame([_junc("chr1", 1000, 1100, "+")])  # 899 bp to gene
    out, _ = cluster_leak(leak, genes=merge_gene_bodies([gtf]))
    assert len(out) == 1
    assert float(out.iloc[0]["nearest_distance_bp"]) > 200


def test_apply_terminal_extents_uses_coverage_and_clips_gene(tmp_path: Path) -> None:
    df = pd.DataFrame(
        [
            {
                "residual_id": "RSDL.1",
                "locus_coord": "chr1:850-1450:+",
                "chrom": "chr1",
                "start": 850,
                "end": 1450,
                "strand": "+",
                "n_junctions": 1,
                "n_exons": 2,
                "length_nt": 300,
                "exon_structure": "850-999,1101-1450",
                "intron_structure": "1000-1100",
                "nearest_gene_name": "Far",
                "nearest_distance_bp": 12000,
                "control_max": 0,
                "treat_sum": 8,
                "treat_n_detected": 3,
                "n_shared_site": 0,
                "status": "unassembled",
            }
        ]
    )
    extents = pd.DataFrame(
        [
            {
                "residual_id": "RSDL.1",
                "chrom": "chr1",
                "strand": "+",
                "intron_start": 1000,
                "intron_end": 1100,
                "left_start": 700,
                "right_end": 1600,
            }
        ]
    )
    out, _ = apply_terminal_extents(df, extents)
    assert int(out.iloc[0]["start"]) == 700
    assert int(out.iloc[0]["end"]) == 1600
    assert "700-999" in str(out.iloc[0]["exon_structure"])
    assert "1101-1600" in str(out.iloc[0]["exon_structure"])
    assert int(out.iloc[0]["length_nt"]) == 800

    genes = {"1": [(650, 720, "+", "Hit", "protein_coding", "G")]}
    left, right = clip_terminals_from_genes("chr1", 700, 1600, 1000, 1100, genes)
    assert left == 721
    assert right == 1600
    clipped, _ = apply_terminal_extents(df, extents, genes)
    assert int(clipped.iloc[0]["n_exons"]) >= 2
    write_residual_gtf(clipped, tmp_path / "residual_cov.gtf")
    for rec in clipped.to_dict(orient="records"):
        for part in str(rec["exon_structure"]).split(","):
            a, b = part.split("-")
            assert not overlaps_gene_body("chr1", int(a), int(b), genes)


def test_terminal_clip_does_not_fall_back_into_antisense_gene() -> None:
    """R1: gene body ends 1 nt before intron. Clip must not restore the 30 nt stub."""
    df = pd.DataFrame(
        [
            {
                "residual_id": "RSDL.1",
                "locus_coord": "chr1:4970-6030:+",
                "chrom": "chr1",
                "start": 4970,
                "end": 6030,
                "strand": "+",
                "n_junctions": 1,
                "n_exons": 2,
                "length_nt": 60,
                "exon_structure": "4970-4999,6001-6030",
                "intron_structure": "5000-6000",
                "nearest_gene_name": "Col",
                "nearest_distance_bp": 0,
                "control_max": 0,
                "treat_sum": 8,
                "treat_n_detected": 3,
                "n_shared_site": 0,
                "status": "unassembled",
            }
        ]
    )
    extents = pd.DataFrame(
        [
            {
                "residual_id": "RSDL.1",
                "left_start": 4970,
                "right_end": 6030,
            }
        ]
    )
    genes = {"1": [(1000, 4999, "-", "Col", "protein_coding", "G")]}
    out, n_deg = apply_terminal_extents(df, extents, genes)
    assert out.empty
    assert n_deg == 1
    empty_ext = pd.DataFrame()
    again, extra = apply_terminal_extents(df, empty_ext, genes)
    assert extra == 0
    assert len(again) == len(df)


def test_residual_gtf_exons_never_overlap_gene_body(tmp_path: Path) -> None:
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(
        'chr1\tX\tgene\t1000\t4999\t.\t-\t.\tgene_id "G"; gene_name "Col"; gene_type "protein_coding";\n',
        encoding="utf-8",
    )
    src = tmp_path / "leak.tsv"
    pd.DataFrame(
        [_junc("chr1", 5000, 6000, "+", nearest_distance_bp=12000, nearest_gene_name="Far")]
    ).to_csv(src, sep="\t", index=False)
    dest = tmp_path / "residual.tsv"
    out = write_residuals(src, dest, gene_gtfs=[gtf])
    genes = merge_gene_bodies([gtf])
    for rec in out.to_dict(orient="records"):
        for part in str(rec["exon_structure"]).split(","):
            if not part:
                continue
            a, b = part.split("-")
            assert not overlaps_gene_body(str(rec["chrom"]), int(a), int(b), genes)


def test_rsdl_ids_are_deterministic() -> None:
    rows = [
        _junc(
            "chr1", 1_000 + 100_000 * i, 1_100 + 100_000 * i, "+", treat_sum=2, treat_n_detected=2
        )
        for i in range(8)
    ]
    leak_a = pd.DataFrame(rows)
    leak_b = pd.DataFrame(list(reversed(rows)))
    a, _ = cluster_leak(leak_a)
    b, _ = cluster_leak(leak_b)
    assert list(a["residual_id"]) == list(b["residual_id"])
    assert list(a["chrom"]) == list(b["chrom"])
    assert list(a["start"]) == list(b["start"])


def test_missing_status_column_is_actionable() -> None:
    from txnova.errors import TxNovaError

    leak = pd.DataFrame([{"chrom": "chr1", "start": 1, "end": 2, "strand": "+"}])
    try:
        select_junctions(leak)
    except TxNovaError as e:
        assert "status" in str(e)
        assert "leak.json" in str(e)
    else:
        raise AssertionError("expected TxNovaError")


def test_knife_runs_after_cluster_not_in_select(tmp_path: Path) -> None:
    gtf = tmp_path / "g.gtf"
    gtf.write_text(
        'chr12\tX\tgene\t1251\t2000\t.\t+\t.\tgene_id "A"; gene_name "Akap"; gene_type "protein_coding";\n',
        encoding="utf-8",
    )
    genes = merge_gene_bodies([gtf])
    leak = pd.DataFrame(
        [_junc("chr12", 1000, 1100, "+", nearest_distance_bp=150, nearest_gene_name="Akap")]
    )
    kept = select_junctions(leak, genes)
    assert len(kept) == 1
    out, _ = cluster_leak(leak, genes=genes)
    assert out.empty


def test_nearest_name_matches_min_distance() -> None:
    leak = pd.DataFrame(
        [
            _junc(
                "chr1",
                1000,
                1100,
                "+",
                nearest_distance_bp=8000,
                nearest_gene_name="Foo",
                treat_sum=3,
            ),
            _junc(
                "chr1",
                1000,
                1300,
                "+",
                nearest_distance_bp=100,
                nearest_gene_name="Bar",
                treat_sum=2,
            ),
        ]
    )
    out, _ = cluster_leak(leak)
    assert len(out) == 1
    assert str(out.iloc[0]["nearest_gene_name"]) == "Bar"
    assert float(out.iloc[0]["nearest_distance_bp"]) == 100


def test_drops_micro_intron_keeps_50nt() -> None:
    leak = pd.DataFrame(
        [
            _junc("chr1", 1000, 1035, "+"),  # 36 nt
            _junc("chr2", 1000, 1049, "+"),  # 50 nt
        ]
    )
    kept = select_junctions(leak)
    assert set(kept["chrom"]) == {"chr2"}


def test_drops_snrna_1bp_keeps_protein_1kb(tmp_path: Path) -> None:
    gtf = tmp_path / "m39.gtf"
    gtf.write_text(
        'chr3\tX\tgene\t1102\t1200\t.\t-\t.\tgene_id "S"; gene_name "Gm24382"; gene_type "snRNA";\n'
        'chr4\tX\tgene\t3000\t4000\t.\t+\t.\tgene_id "P"; gene_name "Far"; gene_type "protein_coding";\n',
        encoding="utf-8",
    )
    genes = merge_gene_bodies([gtf])
    assert close_to_small_rna("chr3", 1000, 1100, genes)
    leak = pd.DataFrame(
        [
            _junc("chr3", 1000, 1100, "+"),  # 1 bp from snRNA
            _junc("chr4", 1000, 1100, "+"),  # 1899 bp from protein
        ]
    )
    out, _ = cluster_leak(leak, genes=genes)
    assert set(out["chrom"]) == {"chr4"}


def test_coverage_stub_terminal_keeps_intron_chain() -> None:
    """A short terminal is missing coverage, not a reason to drop the intron."""
    df = pd.DataFrame(
        [
            {
                "residual_id": "RSDL.1",
                "locus_coord": "chr1:970-1230:+",
                "chrom": "chr1",
                "start": 970,
                "end": 1230,
                "strand": "+",
                "n_junctions": 1,
                "n_exons": 2,
                "length_nt": 60,
                "exon_structure": "970-999,1101-1230",
                "intron_structure": "1000-1100",
                "nearest_gene_name": "Far",
                "nearest_distance_bp": 12000,
                "control_max": 0,
                "control_n_detected": 0,
                "cohort": "silent",
                "treat_sum": 8,
                "treat_n_detected": 3,
                "n_shared_site": 0,
                "status": "unassembled",
            }
        ]
    )
    extents = pd.DataFrame(
        [
            {
                "residual_id": "RSDL.1",
                "left_start": 970,
                "right_end": 1130,
            }
        ]
    )
    out, n_deg = apply_terminal_extents(df, extents)
    assert n_deg == 0
    assert len(out) == 1
    assert str(out.iloc[0]["residual_id"]) == "RSDL.1"
    exons = str(out.iloc[0]["exon_structure"]).split(",")
    assert len(exons) == 2


def test_naming_gtf_union_rejects_known_gene(tmp_path: Path) -> None:
    """Harvest knives use annotation ∪ naming_annotation.

    A mask-and-recall run that leaves the hidden gene in naming_annotation
    will drop every target. Both GTFs must omit the same gene_ids.
    """
    run_gtf = tmp_path / "run.gtf"
    run_gtf.write_text(
        'chr2\tX\tgene\t1\t10\t.\t+\t.\tgene_id "OTHER"; gene_name "Far"; '
        'gene_type "protein_coding";\n',
        encoding="utf-8",
    )
    naming = tmp_path / "naming.gtf"
    naming.write_text(
        'chr1\tX\tgene\t900\t1300\t.\t+\t.\tgene_id "MASK"; '
        'gene_name "Col9a1"; gene_type "protein_coding";\n',
        encoding="utf-8",
    )
    src = tmp_path / "leak.tsv"
    pd.DataFrame([_junc("chr1", 1000, 1100, "+")]).to_csv(src, sep="\t", index=False)
    kept = write_residuals(src, tmp_path / "keep.tsv", gene_gtfs=[run_gtf])
    assert len(kept) == 1
    dropped = write_residuals(src, tmp_path / "drop.tsv", gene_gtfs=[run_gtf, naming])
    assert dropped.empty


def test_missing_or_empty_gtf_is_an_error(tmp_path: Path) -> None:
    missing = tmp_path / "nope.gtf"
    try:
        merge_gene_bodies([missing])
    except TxNovaError as e:
        assert "not found" in str(e)
    else:
        raise AssertionError("expected missing GTF to raise")
    empty = tmp_path / "empty.gtf"
    empty.write_text("# only a comment\n", encoding="utf-8")
    try:
        merge_gene_bodies([empty])
    except TxNovaError as e:
        assert "no gene or transcript" in str(e)
    else:
        raise AssertionError("expected empty GTF to raise")


def test_terminal_clipped_to_contig() -> None:
    chain = [{"chrom": "chr1", "start": 90, "end": 95, "strand": "+"}]
    struct, n_exons, start, end = exon_structure(chain, flank=30, contig_end=100)
    assert end == 100
    assert n_exons == 2
    assert clip_to_contig("chr1", 500, {"chr1": 100}) == 100
    assert clip_to_contig("chrM", 20, {"MT": 16}) == 16
