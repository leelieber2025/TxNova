from __future__ import annotations

from pathlib import Path

import pandas as pd
from txnova.naming import (
    NAMING_ANNOTATED,
    NAMING_NONE,
    NAMING_PROTEIN,
    annotate_table,
    chrom_key,
    classify_locus,
    load_gene_bodies,
)


def test_chrom_key() -> None:
    assert chrom_key("chr8") == "8"
    assert chrom_key("8") == "8"
    assert chrom_key("GL456233.2") == "GL456233.2"
    assert chrom_key("chrM") == "MT"
    assert chrom_key("MT") == "MT"
    assert chrom_key("M") == "MT"


def test_classify_prefers_same_strand_protein_coding(tmp_path: Path) -> None:
    gtf = tmp_path / "g.gtf"
    gtf.write_text(
        'chr4\tHAVANA\tgene\t100\t500\t.\t+\t.\tgene_id "G1"; gene_name "Dhrsx"; gene_type "protein_coding";\n'
        'chr4\tHAVANA\tgene\t200\t400\t.\t-\t.\tgene_id "G2"; gene_name "Other"; gene_type "lncRNA";\n'
        'chr4\tHAVANA\texon\t100\t150\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n',
        encoding="utf-8",
    )
    genes = load_gene_bodies(gtf)
    hit = classify_locus("chr4", 120, 300, "+", genes)
    assert hit["named_overlap"] == NAMING_PROTEIN
    assert hit["named_gene_name"] == "Dhrsx"
    miss = classify_locus("chr4", 900, 950, "+", genes)
    assert miss["named_overlap"] == NAMING_NONE
    lnc = classify_locus("chr8", 1, 10, "+", genes)
    assert lnc["named_overlap"] == NAMING_NONE
    only_lnc = classify_locus("chr4", 200, 400, "-", genes)
    assert only_lnc["named_overlap"] == NAMING_ANNOTATED
    assert only_lnc["named_gene_name"] == "Other"


def test_load_gene_bodies_from_transcript_only(tmp_path: Path) -> None:
    gtf = tmp_path / "refseq.gtf"
    gtf.write_text(
        'chr1\trefGene\ttranscript\t100\t500\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; '
        'gene_name "Col"; gene_biotype "protein_coding";\n'
        'chr1\trefGene\texon\t100\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
        'chr1\trefGene\texon\t400\t500\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n',
        encoding="utf-8",
    )
    genes = load_gene_bodies(gtf)
    recs = genes["1"]
    assert len(recs) == 1
    assert recs[0][0] == 100
    assert recs[0][1] == 500
    assert recs[0][3] == "Col"


def test_annotate_table_adds_columns() -> None:
    genes = {
        "1": [(1000, 2000, "+", "Ccsap", "protein_coding", "ENS1")],
    }
    df = pd.DataFrame(
        [
            {
                "locus_id": "L1",
                "chrom": "chr1",
                "start": 1500,
                "end": 1800,
                "strand": "+",
            },
            {
                "locus_id": "L2",
                "chrom": "chr1",
                "start": 5000,
                "end": 5100,
                "strand": "+",
            },
        ]
    )
    out = annotate_table(df, genes)
    assert list(out.loc[out.locus_id == "L1", "named_gene_name"]) == ["Ccsap"]
    assert list(out.loc[out.locus_id == "L2", "named_overlap"]) == [NAMING_NONE]
