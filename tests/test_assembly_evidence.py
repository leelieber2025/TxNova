from __future__ import annotations

from pathlib import Path

import pandas as pd
from txnova.assembly_evidence import attach_assembly_evidence
from txnova.samples import SampleRow


def test_assembly_evidence_treat_only(tmp_path: Path) -> None:
    asm = tmp_path / "assembly"
    asm.mkdir()
    (asm / "c1.gtf").write_text(
        'chr1\tStringTie\ttranscript\t10\t80\t.\t+\t.\tgene_id "STRG.1"; transcript_id "STRG.1.1"; cov "3.0";\n',
        encoding="utf-8",
    )
    (asm / "t1.gtf").write_text(
        'chr1\tStringTie\ttranscript\t1000\t1200\t.\t+\t.\tgene_id "STRG.9"; transcript_id "STRG.9.1"; cov "8.8";\n'
        'chr1\tStringTie\ttranscript\t1000\t1200\t.\t+\t.\tgene_id "STRG.9"; transcript_id "STRG.9.2"; cov "2.1";\n',
        encoding="utf-8",
    )
    (asm / "t2.gtf").write_text(
        'chr1\tStringTie\ttranscript\t980\t1250\t.\t+\t.\tgene_id "STRG.8"; transcript_id "STRG.8.1"; cov "12.0";\n',
        encoding="utf-8",
    )
    cand = tmp_path / "candidates.tsv"
    cand.write_text(
        "locus_id\tchrom\tstart\tend\tstrand\nL1\tchr1\t1000\t1200\t+\n",
        encoding="utf-8",
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
    attach_assembly_evidence(cand, rows, asm)
    df = pd.read_csv(cand, sep="\t")
    assert df.iloc[0]["assembled_in"] == "t1,t2"
    assert int(df.iloc[0]["n_assembled_control"]) == 0
    assert int(df.iloc[0]["n_assembled_treat"]) == 2
    assert bool(df.iloc[0]["t1_assembled"])
    assert not bool(df.iloc[0]["c1_assembled"])
    assert int(df.iloc[0]["t1_asm_n_isoforms"]) == 2
    assert float(df.iloc[0]["t1_asm_max_cov"]) == 8.8
    assert df.iloc[0]["t2_asm_span"] == "980-1250"
    assert int(df.iloc[0]["asm_span_max_delta_bp"]) == 50
