"""TxNova assembler: annotation + residual splice loci.

StringTie 3.0.3 (gpertea/stringtie) does three expensive things we do not
need for treat-specific intergenic discovery:

- bundle.cpp / BundleData: pack every overlapping read on a locus
- rlink.cpp: per-sample splice graph + flow (isoform abundance)
- tmerge.cpp: --merge consensus → hundreds of thousands of MSTRG rows

Taken from StringTie, implemented here / in residual.py:

- CJunction: intron is the node (leak CIGAR N, treat-recurrent)
- coverage trim of terminals (StringTie -t default; our treat-depth walk)
- guide annotation as the known universe (StringTie -G), not as a rebuild

assemble() does not scan the BAM. merge() writes gene/transcript/exon
rows from the reference GTF (CDS/UTR dropped; rust parse ignores them
anyway). Novel models are residual.gtf, streamed into universe.gtf.
"""

from __future__ import annotations

from pathlib import Path

from txnova.errors import TxNovaError
from txnova.logging import get_logger

log = get_logger("txnova.assemble")

ASSEMBLER_NAME = "txnova-residual"
ASSEMBLER_VERSION = "txnova-residual-2"

# rust parse_gtf only keeps these; CDS/UTR/codon rows are dead weight.
MODEL_FEATURES = frozenset({"gene", "transcript", "exon"})


class Assembler:
    """Only contact surface for assembly. No StringTie subprocess."""

    def __init__(self, binary: str | None = None) -> None:
        self.binary_name = binary or ASSEMBLER_NAME

    def resolve_binary(self, name: str | None = None) -> Path:
        return Path(name or self.binary_name)

    def version(self) -> str:
        return ASSEMBLER_VERSION

    def assemble(
        self,
        bam: Path,
        gtf_ref: Path,
        out_gtf: Path,
        *,
        strandedness: str,
        threads: int,
        extra_args: list[str],
    ) -> None:
        # Per-sample transcriptome reconstruction is StringTie's cost and
        # not our discovery engine. Residual leak is multi-sample.
        del bam, gtf_ref, strandedness, threads, extra_args
        out_gtf.parent.mkdir(parents=True, exist_ok=True)
        out_gtf.write_text(
            f"# {ASSEMBLER_NAME}: no per-sample reconstruct (see residual.gtf)\n",
            encoding="utf-8",
        )
        out_gtf.with_suffix(".log").write_text(
            "skipped StringTie-style per-sample assemble\n", encoding="utf-8"
        )

    def merge(
        self,
        sample_gtfs: list[Path],
        gtf_ref: Path,
        out_gtf: Path,
        *,
        extra_args: list[str],
    ) -> None:
        del sample_gtfs, extra_args
        if not gtf_ref.is_file():
            raise TxNovaError(f"annotation GTF missing: {gtf_ref}")
        out_gtf.parent.mkdir(parents=True, exist_ok=True)
        n_tx = _write_model_gtf(gtf_ref, out_gtf)
        if n_tx == 0:
            raise TxNovaError(f"annotation GTF has 0 transcripts: {gtf_ref}")
        out_gtf.with_suffix(".log").write_text(
            f"slim gene/transcript/exon {gtf_ref} → {out_gtf} ({n_tx} transcripts)\n",
            encoding="utf-8",
        )
        log.info("merge: slim annotation %s → %s (%s transcripts)", gtf_ref, out_gtf, n_tx)


def is_model_gtf_line(line: str) -> bool:
    """Keep comments and gene/transcript/exon. Drop CDS/UTR/codon."""
    if not line.strip() or line.startswith("#"):
        return True
    cols = line.split("\t", 3)
    return len(cols) >= 3 and cols[2] in MODEL_FEATURES


def _write_model_gtf(src: Path, dest: Path) -> int:
    """Stream gene/transcript/exon into a new file. Never share the source inode."""
    if dest.exists() or dest.is_symlink():
        dest.unlink()
    n_tx = 0
    with (
        src.open(encoding="utf-8", errors="replace") as fh,
        dest.open("w", encoding="utf-8") as out,
    ):
        for line in fh:
            if not is_model_gtf_line(line):
                continue
            out.write(line)
            cols = line.split("\t", 3)
            if len(cols) >= 3 and cols[2] == "transcript":
                n_tx += 1
    return n_tx


def _n_transcripts(gtf: Path) -> int:
    if not gtf.is_file():
        return 0
    n = 0
    with gtf.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.split("\t", 3)
            if len(cols) >= 3 and cols[2] == "transcript":
                n += 1
    return n
