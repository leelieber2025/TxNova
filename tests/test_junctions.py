from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from txnova import _core

FIXTURES = Path(__file__).resolve().parent / "fixtures"


def test_structure_scan_junction_and_bridge(tmp_path: Path) -> None:
    bam = FIXTURES / "junc_bridge.bam"
    merged = tmp_path / "merged.gtf"
    merged.write_text(
        'chr1\tX\ttranscript\t400\t550\t.\t+\t.\tgene_id "MSTRG.1"; transcript_id "MSTRG.1.1";\n'
        'chr1\tX\texon\t400\t450\t.\t+\t.\tgene_id "MSTRG.1"; transcript_id "MSTRG.1.1";\n'
        'chr1\tX\texon\t500\t550\t.\t+\t.\tgene_id "MSTRG.1"; transcript_id "MSTRG.1.1";\n',
        encoding="utf-8",
    )
    ref = tmp_path / "ref.gtf"
    ref.write_text(
        'chr1\tX\tgene\t10\t250\t.\t+\t.\tgene_id "G1"; gene_name "Car";\n'
        'chr1\tX\ttranscript\t10\t250\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "Car";\n'
        'chr1\tX\texon\t10\t80\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
        'chr1\tX\texon\t200\t250\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n',
        encoding="utf-8",
    )
    reps = tmp_path / "reps.tsv"
    reps.write_text("locus_id\ttranscript_id\nMSTRG.1\tMSTRG.1.1\n", encoding="utf-8")
    class_tsv = tmp_path / "class.tsv"
    class_tsv.write_text("transcript_id\tgene_id\tclass\nMSTRG.1.1\tMSTRG.1\tu\n")
    out = tmp_path / "struct.tsv"
    payload = json.dumps(
        {
            "samples": [
                {"sample_id": "t1", "bam": str(bam), "dup_flag_seen": False},
            ]
        }
    )
    cfg = json.dumps(
        {
            "strandedness": "unstranded",
            "min_mapq": 10,
            "skip_duplicate": "never",
            "discontinuity_window_bp": 50,
            "discontinuity_valley_bp": 50,
            "threads": 1,
        }
    )
    _core.structure_scan(
        str(merged),
        str(ref),
        str(class_tsv),
        str(reps),
        str(FIXTURES / "genome.fa"),
        payload,
        str(out),
        cfg,
    )
    df = pd.read_csv(out, sep="\t")
    assert len(df) == 1
    assert str(df.iloc[0]["junction_support"]) == "1"
    assert int(df.iloc[0]["bridge_read_count"]) == 1
    assert str(df.iloc[0]["nearest_gene_name"]) == "Car"
