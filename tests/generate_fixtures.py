#!/usr/bin/env python3
"""Generate committed test fixtures. Not run in CI."""

from __future__ import annotations

from pathlib import Path

import pysam

ROOT = Path(__file__).resolve().parent / "fixtures"
SEQ = "ACGT" * 250  # 1000 bp
CHROM = "chr1"
LN = len(SEQ)


def write_fasta(path: Path, name: str, seq: str) -> None:
    path.write_text(f">{name}\n{seq}\n", encoding="ascii")
    # Manual .fai: name length offset linebases linewidth
    # >name\n then seq\n  → offset = 1 + len(name) + 1
    offset = 1 + len(name) + 1
    fai = path.with_name(path.name + ".fai")
    fai.write_text(f"{name}\t{len(seq)}\t{offset}\t{len(seq)}\t{len(seq) + 1}\n", encoding="ascii")


def write_gtf(path: Path, seqname: str) -> None:
    # one gene, two exons — seqname only matters for preflight
    path.write_text(
        f'{seqname}\tHAVANA\tgene\t10\t200\t.\t+\t.\tgene_id "ENSG1"; gene_name "G1";\n'
        f'{seqname}\tHAVANA\ttranscript\t10\t200\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1";\n'
        f'{seqname}\tHAVANA\texon\t10\t80\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1";\n'
        f'{seqname}\tHAVANA\texon\t120\t200\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1";\n',
        encoding="ascii",
    )


def make_header(pg: list[dict]) -> dict:
    return {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": CHROM, "LN": LN}],
        "PG": pg,
    }


def write_se_bam(path: Path, pg: list[dict], n_reads: int = 4) -> None:
    header = make_header(pg)
    path.parent.mkdir(parents=True, exist_ok=True)
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for i in range(n_reads):
            a = pysam.AlignedSegment()
            a.query_name = f"r{i}"
            a.query_sequence = "ACGT" * 10
            a.flag = 0
            a.reference_id = 0
            a.reference_start = 10 + i * 5
            a.mapping_quality = 255
            a.cigar = ((0, 40),)
            a.query_qualities = pysam.qualitystring_to_array("I" * 40)
            out.write(a)
    pysam.index(str(path))


def main() -> None:
    ROOT.mkdir(parents=True, exist_ok=True)
    write_fasta(ROOT / "genome.fa", CHROM, SEQ)
    write_gtf(ROOT / "genes.gtf", CHROM)
    write_gtf(ROOT / "genes_nochr.gtf", "1")

    star_samtools = [
        {"ID": "STAR", "PN": "STAR", "VN": "2.7.11a"},
        {"ID": "samtools", "PN": "samtools", "PP": "STAR", "VN": "1.20"},
    ]
    bowtie = [{"ID": "bowtie2", "PN": "bowtie2", "VN": "2.5.3"}]
    star_bowtie = [
        {"ID": "STAR", "PN": "STAR", "VN": "2.7.11a"},
        {"ID": "bowtie2", "PN": "bowtie2"},
    ]

    for name in ("ctrl_1", "ctrl_2", "treat_1", "treat_2"):
        write_se_bam(ROOT / f"{name}.bam", star_samtools)

    write_se_bam(ROOT / "bowtie2.bam", bowtie)
    write_se_bam(ROOT / "star_bowtie2.bam", star_bowtie)

    src = ROOT / "ctrl_1.bam"
    truncated = ROOT / "truncated.bam"
    data = src.read_bytes()
    truncated.write_bytes(data[:-40])
    # no index on purpose

    noindex = ROOT / "noindex.bam"
    noindex.write_bytes(src.read_bytes())
    # deliberately no .bai

    ok_sheet = ROOT / "samples_ok.tsv"
    ok_sheet.write_text(
        "sample_id\tbam\tgroup\tstrandedness\treplicate\n"
        + "".join(
            f"{s}\t{s}.bam\t{'control' if s.startswith('ctrl') else 'treat'}\trf\t{1 if s.endswith('1') else 2}\n"
            for s in ("ctrl_1", "ctrl_2", "treat_1", "treat_2")
        ),
        encoding="utf-8",
    )

    (ROOT / "config_ok.yaml").write_text(
        """\
species: mouse
output_dir: ./txnova_out
threads: 0
genome:
  fasta: genome.fa
  annotation: genes.gtf
  annotation_source: GENCODE
  annotation_version: M39
  assembly: GRCm39
samples: samples_ok.tsv
""",
        encoding="utf-8",
    )
    _write_quant_and_junc_bams()
    print(f"wrote fixtures under {ROOT}")


def _aln(
    name: str,
    start0: int,
    length: int,
    flag: int = 0,
    reverse: bool = False,
    cigar=None,
    nh: int | None = 1,
) -> pysam.AlignedSegment:
    a = pysam.AlignedSegment()
    a.query_name = name
    a.query_sequence = "A" * length
    a.flag = flag | (16 if reverse else 0)
    a.reference_id = 0
    a.reference_start = start0
    a.mapping_quality = 255
    a.cigar = cigar if cigar is not None else ((0, length),)
    a.query_qualities = pysam.qualitystring_to_array("I" * length)
    if nh is not None:
        a.set_tag("NH", nh)
    return a


def _write_bam(path: Path, records: list[pysam.AlignedSegment]) -> None:
    header = make_header([{"ID": "STAR", "PN": "STAR", "VN": "2.7.11a"}])
    path.parent.mkdir(parents=True, exist_ok=True)
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for rec in records:
            out.write(rec)
    pysam.index(str(path))


def _write_quant_and_junc_bams() -> None:
    _write_bam(ROOT / "quant_intron_only.bam", [_aln("r0", 200, 40, nh=None)])
    _write_bam(ROOT / "quant_exon_hit.bam", [_aln("r0", 120, 30, nh=None)])
    _write_bam(ROOT / "quant_cross_locus.bam", [_aln("r0", 160, 20, nh=None)])
    _write_bam(ROOT / "quant_tpm_novel.bam", [_aln("r0", 520, 20, nh=None)])
    _write_bam(
        ROOT / "quant_tpm_both.bam",
        [_aln("r0", 120, 20, nh=None), _aln("r1", 520, 20, nh=None)],
    )
    _write_bam(ROOT / "quant_secondary.bam", [_aln("r0", 120, 30, flag=256, nh=None)])

    bridge = _aln("bridge1", 199, 102, cigar=((0, 51), (3, 149), (0, 51)))
    mate = _aln("bridge1", 199, 102, flag=128, cigar=((0, 51), (3, 149), (0, 51)))
    junc = _aln("junc1", 399, 101, cigar=((0, 51), (3, 49), (0, 50)))
    _write_bam(ROOT / "junc_bridge.bam", [bridge, mate, junc])


if __name__ == "__main__":
    main()
