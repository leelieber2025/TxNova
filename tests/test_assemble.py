from __future__ import annotations

from pathlib import Path

from txnova.assemble import ASSEMBLER_VERSION, Assembler, _n_transcripts


def test_assemble_does_not_call_stringtie(tmp_path: Path) -> None:
    dest = tmp_path / "s1.gtf"
    Assembler().assemble(
        tmp_path / "missing.bam",
        tmp_path / "ref.gtf",
        dest,
        strandedness="rf",
        threads=8,
        extra_args=["--rf"],
    )
    text = dest.read_text()
    assert text.startswith("# txnova-residual")
    assert _n_transcripts(dest) == 0
    assert dest.with_suffix(".log").is_file()


def test_merge_writes_slim_model_gtf(tmp_path: Path) -> None:
    ref = tmp_path / "ref.gtf"
    ref.write_text(
        'chr1\tHAVANA\tgene\t1\t10\t.\t+\t.\tgene_id "G";\n'
        'chr1\tHAVANA\ttranscript\t1\t10\t.\t+\t.\tgene_id "G"; transcript_id "T";\n'
        'chr1\tHAVANA\texon\t1\t10\t.\t+\t.\tgene_id "G"; transcript_id "T";\n'
        'chr1\tHAVANA\tCDS\t3\t8\t.\t+\t.\tgene_id "G"; transcript_id "T";\n',
        encoding="utf-8",
    )
    dest = tmp_path / "merged.gtf"
    Assembler().merge([], ref, dest, extra_args=[])
    text = dest.read_text()
    assert "\ttranscript\t" in text
    assert "\texon\t" in text
    assert "\tCDS\t" not in text
    assert not dest.is_symlink()
    assert dest.stat().st_ino != ref.stat().st_ino
    assert _n_transcripts(dest) == 1
    assert Assembler().version() == ASSEMBLER_VERSION
