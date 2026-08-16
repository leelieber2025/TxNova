"""Residual splice loci: cluster leak junctions into exon models.

Borrows LeafCutter's connected-component idea (shared splice sites) and
adds constitutive chaining (plausible exon between adjacent introns).
This is the novel-model half of the assembler. Does not call
LeafCutter/Portcullis binaries.
"""

from __future__ import annotations

import json
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import pandas as pd

from txnova.assemble import is_model_gtf_line
from txnova.errors import TxNovaError
from txnova.logging import get_logger
from txnova.naming import chrom_key, load_gene_bodies
from txnova.orphan import parse_exons

log = get_logger("txnova.residual")

# gene body: (start, end, strand, name, type, id)
GeneBodies = dict[str, list[tuple[int, int, str, str, str, str]]]

RESIDUAL_FILENAME = "residual.tsv"
RESIDUAL_GTF = "residual.gtf"
RESIDUAL_STATS = "residual.stats.json"
UNIVERSE_GTF = "universe.gtf"
RESIDUAL_COLUMNS = [
    "residual_id",
    "locus_coord",
    "chrom",
    "start",
    "end",
    "strand",
    "n_junctions",
    "n_exons",
    "length_nt",
    "exon_structure",
    "intron_structure",
    "nearest_gene_name",
    "nearest_distance_bp",
    "control_max",
    "treat_sum",
    "treat_n_detected",
    "n_shared_site",
    "status",
]

# Internal exon between two residual introns (not first/last exon).
MIN_EXON_NT = 30
MAX_EXON_NT = 20_000
# Fallback stub if there is no treat BAM / no coverage at the splice.
MIN_TERMINAL_NT = 30
TERMINAL_FLANK_NT = MIN_TERMINAL_NT
# Treat-coverage walk from each splice site (Cufflinks / StringTie style).
MAX_TERMINAL_NT = 2_000
MIN_TERMINAL_DEPTH = 1
MAX_TERMINAL_GAP_NT = 20
# Same-strand, comprehensive GTF: drop only if stuck to a gene (extension).
# Not a 5 kb gate. Opposite strand is ignored here.
CLOSE_SAME_STRAND_NT = 200


def _true(v: Any) -> bool:
    return str(v).strip().lower() in {"true", "1", "yes"}


def _num(v: Any) -> float | None:
    if v is None or v == "" or (isinstance(v, float) and pd.isna(v)):
        return None
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def merge_gene_bodies(gtfs: list[Path] | None) -> GeneBodies:
    out: GeneBodies = defaultdict(list)
    for raw in gtfs or []:
        p = Path(raw)
        if not p.is_file():
            continue
        for chrom, recs in load_gene_bodies(p).items():
            out[chrom].extend(recs)
    return dict(out)


def overlaps_gene_body(
    chrom: str,
    start: int,
    end: int,
    genes: GeneBodies,
    *,
    strand: str | None = None,
) -> bool:
    """Either-strand overlap unless `strand` is set (then same-strand only)."""
    for gs, ge, gstrand, *_ in genes.get(chrom_key(chrom), []):
        if strand is not None and gstrand != strand:
            continue
        if start <= ge and gs <= end:
            return True
    return False


def nearest_same_strand(
    chrom: str,
    start: int,
    end: int,
    strand: str,
    genes: GeneBodies,
) -> tuple[str, int | None]:
    best_name = ""
    best_d: int | None = None
    for gs, ge, gstrand, name, *_ in genes.get(chrom_key(chrom), []):
        if gstrand != strand:
            continue
        if start <= ge and gs <= end:
            return name, 0
        if end < gs:
            d = gs - end - 1
        else:
            d = start - ge - 1
        if best_d is None or d < best_d:
            best_d = d
            best_name = name
    return best_name, best_d


def select_junctions(leak: pd.DataFrame, genes: GeneBodies | None = None) -> pd.DataFrame:
    """Unassembled, not overlapping a gene body on either strand.

    The 200 nt same-strand knife is contract R step 7 (after clustering),
    not applied here.
    """
    if leak.empty:
        return leak
    if "status" not in leak.columns:
        raise TxNovaError("leak.tsv missing column status; delete stamps/leak.json and rerun")
    keep = leak["status"].astype(str).eq("unassembled")
    if "overlaps_gene" in leak.columns:
        keep = keep & ~leak["overlaps_gene"].map(_true)
    src = leak.loc[keep].copy()
    if src.empty or not genes:
        return src
    keep_rows = []
    for rec in src.to_dict(orient="records"):
        chrom, start, end = str(rec["chrom"]), int(rec["start"]), int(rec["end"])
        keep_rows.append(not overlaps_gene_body(chrom, start, end, genes))
    return src.iloc[[i for i, ok in enumerate(keep_rows) if ok]].copy()


def same_locus(
    a: dict,
    b: dict,
    *,
    min_exon: int = MIN_EXON_NT,
    max_exon: int = MAX_EXON_NT,
) -> bool:
    if a["chrom"] != b["chrom"] or a["strand"] != b["strand"]:
        return False
    if a["start"] == b["start"] or a["end"] == b["end"]:
        return True
    if a["start"] <= b["end"] and b["start"] <= a["end"]:
        return True
    if a["end"] < b["start"]:
        gap = b["start"] - a["end"] - 1
    else:
        gap = a["start"] - b["end"] - 1
    return min_exon <= gap <= max_exon


class _UF:
    def __init__(self, n: int) -> None:
        self.p = list(range(n))

    def find(self, i: int) -> int:
        while self.p[i] != i:
            self.p[i] = self.p[self.p[i]]
            i = self.p[i]
        return i

    def union(self, i: int, j: int) -> None:
        ri, rj = self.find(i), self.find(j)
        if ri != rj:
            self.p[rj] = ri


def cluster_junctions(
    rows: list[dict],
    *,
    min_exon: int = MIN_EXON_NT,
    max_exon: int = MAX_EXON_NT,
) -> list[list[dict]]:
    if not rows:
        return []
    buckets: dict[tuple[str, str], list[int]] = defaultdict(list)
    for i, r in enumerate(rows):
        buckets[(str(r["chrom"]), str(r["strand"]))].append(i)
    uf = _UF(len(rows))
    for idxs in buckets.values():
        by_start: dict[int, list[int]] = defaultdict(list)
        by_end: dict[int, list[int]] = defaultdict(list)
        for i in idxs:
            by_start[int(rows[i]["start"])].append(i)
            by_end[int(rows[i]["end"])].append(i)
        for group in (*by_start.values(), *by_end.values()):
            head = group[0]
            for j in group[1:]:
                uf.union(head, j)
        ordered = sorted(idxs, key=lambda i: (int(rows[i]["start"]), int(rows[i]["end"])))
        for a_pos, ia in enumerate(ordered):
            a_end = int(rows[ia]["end"])
            limit = a_end + max_exon + 1
            for ib in ordered[a_pos + 1 :]:
                if int(rows[ib]["start"]) > limit:
                    break
                if same_locus(rows[ia], rows[ib], min_exon=min_exon, max_exon=max_exon):
                    uf.union(ia, ib)
    groups: dict[int, list[dict]] = defaultdict(list)
    for i, r in enumerate(rows):
        groups[uf.find(i)].append(r)
    return list(groups.values())


def representative_introns(members: list[dict]) -> list[dict]:
    """Non-overlapping chain; prefer higher treat_sum, then shorter intron."""
    ordered = sorted(
        members,
        key=lambda r: (-float(r.get("treat_sum") or 0), int(r["start"]), int(r["end"])),
    )
    chain: list[dict] = []
    for r in ordered:
        s, e = int(r["start"]), int(r["end"])
        if any(s <= int(c["end"]) and int(c["start"]) <= e for c in chain):
            continue
        chain.append(r)
    chain.sort(key=lambda r: int(r["start"]))
    return chain


def exon_structure(
    chain: list[dict], *, flank: int = TERMINAL_FLANK_NT
) -> tuple[str, int, int, int]:
    """Return exon_structure, n_exons, locus_start, locus_end (1-based inclusive)."""
    if not chain:
        return "", 0, 0, 0
    first = int(chain[0]["start"])
    last = int(chain[-1]["end"])
    start = max(1, first - flank)
    exons: list[tuple[int, int]] = [(start, first - 1)]
    for a, b in zip(chain, chain[1:]):
        exons.append((int(a["end"]) + 1, int(b["start"]) - 1))
    exons.append((last + 1, last + flank))
    exons = [(s, e) for s, e in exons if e >= s]
    struct = ",".join(f"{s}-{e}" for s, e in exons)
    return struct, len(exons), start, last + flank


def _rebuild_exons(intron_structure: str, left: int, right: int) -> tuple[str, int, int, int, int]:
    """Spliced model from intron chain + terminal bounds. length last."""
    introns = parse_exons(intron_structure)
    if not introns:
        return "", 0, 0, 0, 0
    first, last = introns[0][0], introns[-1][1]
    exons: list[tuple[int, int]] = []
    if left < first:
        exons.append((max(1, left), first - 1))
    for a, b in zip(introns, introns[1:]):
        exons.append((a[1] + 1, b[0] - 1))
    if right > last:
        exons.append((last + 1, right))
    exons = [(s, e) for s, e in exons if e >= s]
    if not exons:
        return "", 0, first, last, 0
    loc_s, loc_e = exons[0][0], exons[-1][1]
    struct = ",".join(f"{s}-{e}" for s, e in exons)
    length = sum(e - s + 1 for s, e in exons)
    return struct, len(exons), loc_s, loc_e, length


def clip_terminals_from_genes(
    chrom: str,
    left: int,
    right: int,
    intron_s: int,
    intron_e: int,
    genes: GeneBodies | None,
) -> tuple[int, int]:
    """Do not let a terminal walk into any gene body (either strand)."""
    if not genes:
        return left, right
    left_end = intron_s - 1
    right_start = intron_e + 1
    for gs, ge, *_ in genes.get(chrom_key(chrom), []):
        if left <= ge and gs <= left_end:
            left = max(left, ge + 1)
        if right_start <= ge and gs <= right:
            right = min(right, gs - 1)
    return left, right


def apply_terminal_extents(
    df: pd.DataFrame,
    extents: pd.DataFrame,
    genes: GeneBodies | None = None,
    *,
    min_terminal: int = MIN_TERMINAL_NT,
) -> tuple[pd.DataFrame, int]:
    """Replace stub terminals with coverage walks, clipped off gene bodies."""
    if df.empty or extents is None or extents.empty:
        return df, 0
    by_id = extents.set_index("residual_id")
    rows: list[dict] = []
    n_degenerate = 0
    for rec in df.to_dict(orient="records"):
        rid = str(rec["residual_id"])
        introns = parse_exons(rec.get("intron_structure"))
        if not introns:
            rows.append(rec)
            continue
        intron_s, intron_e = introns[0][0], introns[-1][1]
        stub_left = max(1, intron_s - min_terminal)
        stub_right = intron_e + min_terminal
        if rid in by_id.index:
            hit = by_id.loc[rid]
            raw_l = pd.to_numeric(hit["left_start"], errors="coerce")
            raw_r = pd.to_numeric(hit["right_end"], errors="coerce")
            left = stub_left if pd.isna(raw_l) else int(raw_l)
            right = stub_right if pd.isna(raw_r) else int(raw_r)
        else:
            left, right = stub_left, stub_right
        left = min(left, stub_left)
        right = max(right, stub_right)
        left, right = clip_terminals_from_genes(
            str(rec["chrom"]), left, right, intron_s, intron_e, genes
        )
        # No room after clip: drop that terminal. Never fall back into a gene.
        if left > intron_s - 1:
            left = intron_s
        if right < intron_e + 1:
            right = intron_e
        struct, n_exons, loc_s, loc_e, length = _rebuild_exons(
            str(rec["intron_structure"]), left, right
        )
        if n_exons < 2:
            n_degenerate += 1
            continue
        rec["start"] = loc_s
        rec["end"] = loc_e
        rec["n_exons"] = n_exons
        rec["length_nt"] = length
        rec["exon_structure"] = struct
        rec["locus_coord"] = f"{rec['chrom']}:{loc_s}-{loc_e}:{rec['strand']}"
        rows.append(rec)
    if not rows:
        return pd.DataFrame(columns=RESIDUAL_COLUMNS), n_degenerate
    return pd.DataFrame(rows)[RESIDUAL_COLUMNS], n_degenerate


def _cluster_record(members: list[dict], residual_id: str) -> dict:
    chain = representative_introns(members)
    struct, n_exons, loc_s, loc_e = exon_structure(chain)
    nearest_name, nearest_d = _nearest_from_members(members)
    treat_n = max(int(m.get("treat_n_detected") or 0) for m in members)
    treat_sum = int(sum(int(m.get("treat_sum") or 0) for m in members))
    control_max = max(int(m.get("control_max") or 0) for m in members)
    return {
        "residual_id": residual_id,
        "chrom": members[0]["chrom"],
        "start": loc_s,
        "end": loc_e,
        "strand": members[0]["strand"],
        "n_junctions": len(members),
        "n_exons": n_exons,
        "length_nt": sum(e - s + 1 for s, e in (map(int, p.split("-")) for p in struct.split(",")))
        if struct
        else 0,
        "exon_structure": struct,
        "intron_structure": ",".join(f"{int(c['start'])}-{int(c['end'])}" for c in chain),
        "nearest_gene_name": nearest_name,
        "nearest_distance_bp": nearest_d,
        "control_max": control_max,
        "treat_sum": treat_sum,
        "treat_n_detected": treat_n,
        "n_shared_site": _n_shared_site(members),
        "status": "unassembled",
    }


def _nearest_from_members(members: list[dict]) -> tuple[str, object]:
    """Name and distance from the same intron (argmin of nearest_distance_bp)."""
    best_d: float | None = None
    best_name = ""
    for m in members:
        d = _num(m.get("nearest_distance_bp"))
        if d is None:
            continue
        name = str(m.get("nearest_gene_name") or "")
        if name in {"", "nan", "NA"}:
            name = ""
        if best_d is None or d < best_d:
            best_d = d
            best_name = name
    return best_name, (best_d if best_d is not None else "")


def _n_shared_site(members: list[dict]) -> int:
    """Count splice sites used by more than one intron (not pairwise)."""
    donors = Counter(int(m["start"]) for m in members)
    accs = Counter(int(m["end"]) for m in members)
    return sum(1 for c in donors.values() if c > 1) + sum(1 for c in accs.values() if c > 1)


def _intron_span(members: list[dict]) -> tuple[int, int]:
    chain = representative_introns(members)
    if not chain:
        return 0, 0
    return min(int(c["start"]) for c in chain), max(int(c["end"]) for c in chain)


def cluster_leak(
    leak: pd.DataFrame,
    *,
    min_exon: int = MIN_EXON_NT,
    max_exon: int = MAX_EXON_NT,
    genes: GeneBodies | None = None,
) -> tuple[pd.DataFrame, int]:
    src = select_junctions(leak, genes)
    if src.empty:
        return pd.DataFrame(columns=RESIDUAL_COLUMNS), 0
    rows = src.to_dict(orient="records")
    clusters = cluster_junctions(rows, min_exon=min_exon, max_exon=max_exon)
    out_rows: list[dict] = []
    n_degenerate = 0
    for members in clusters:
        rec = _cluster_record(members, "")
        # Intron span only. Terminal flanks must not kill a nearby locus.
        if genes:
            span_s, span_e = _intron_span(members)
            chrom, strand = str(rec["chrom"]), str(rec["strand"])
            if overlaps_gene_body(chrom, span_s, span_e, genes):
                continue
            _name, close_d = nearest_same_strand(chrom, span_s, span_e, strand, genes)
            if close_d is not None and close_d < CLOSE_SAME_STRAND_NT:
                continue
            clipped = _clip_record_terminals(rec, genes)
            if clipped is None:
                n_degenerate += 1
                continue
            rec = clipped
        rec["locus_coord"] = f"{rec['chrom']}:{rec['start']}-{rec['end']}:{rec['strand']}"
        out_rows.append(rec)
    if not out_rows:
        return pd.DataFrame(columns=RESIDUAL_COLUMNS), n_degenerate
    out = pd.DataFrame(out_rows)
    out = out.sort_values(
        ["n_junctions", "treat_sum", "treat_n_detected", "chrom", "start"],
        ascending=[False, False, False, True, True],
        kind="mergesort",
    )
    out["residual_id"] = [f"RSDL.{i}" for i in range(1, len(out) + 1)]
    return out[RESIDUAL_COLUMNS], n_degenerate


def _clip_record_terminals(rec: dict, genes: GeneBodies) -> dict | None:
    """Crop stub/coverage terminals off gene bodies. Degenerate if no room."""
    introns = parse_exons(rec.get("intron_structure"))
    exons = parse_exons(rec.get("exon_structure"))
    if not introns or not exons:
        return rec
    intron_s, intron_e = introns[0][0], introns[-1][1]
    left = exons[0][0] if exons[0][1] == intron_s - 1 else intron_s
    right = exons[-1][1] if exons[-1][0] == intron_e + 1 else intron_e
    left, right = clip_terminals_from_genes(
        str(rec["chrom"]), left, right, intron_s, intron_e, genes
    )
    if left > intron_s - 1:
        left = intron_s
    if right < intron_e + 1:
        right = intron_e
    struct, n_exons, loc_s, loc_e, length = _rebuild_exons(
        str(rec["intron_structure"]), left, right
    )
    if n_exons < 2:
        return None
    rec["start"] = loc_s
    rec["end"] = loc_e
    rec["n_exons"] = n_exons
    rec["length_nt"] = length
    rec["exon_structure"] = struct
    rec["locus_coord"] = f"{rec['chrom']}:{loc_s}-{loc_e}:{rec['strand']}"
    return rec


def _write_empty(dest: Path) -> pd.DataFrame:
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text("\t".join(RESIDUAL_COLUMNS) + "\n", encoding="utf-8")
    gtf = dest.with_name(RESIDUAL_GTF)
    gtf.write_text("", encoding="utf-8")
    _write_residual_stats(dest, 0, 0)
    return pd.DataFrame(columns=RESIDUAL_COLUMNS)


def write_residual_gtf(df: pd.DataFrame, dest: Path) -> None:
    """StringTie-shaped GTF: gene_id=RSDL.n, transcript_id=RSDL.n.1."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if df is None or df.empty:
        dest.write_text("", encoding="utf-8")
        return
    lines: list[str] = []
    for rec in df.to_dict(orient="records"):
        rid = str(rec["residual_id"])
        chrom = str(rec["chrom"])
        strand = str(rec["strand"])
        exons = parse_exons(rec.get("exon_structure"))
        if not exons:
            continue
        exons = sorted(exons)
        gs, ge = exons[0][0], exons[-1][1]
        tid = f"{rid}.1"
        gattr = f'gene_id "{rid}";'
        tattr = f'gene_id "{rid}"; transcript_id "{tid}";'
        lines.append(f"{chrom}\ttxnova\tgene\t{gs}\t{ge}\t.\t{strand}\t.\t{gattr}")
        lines.append(f"{chrom}\ttxnova\ttranscript\t{gs}\t{ge}\t.\t{strand}\t.\t{tattr}")
        for s, e in exons:
            lines.append(f"{chrom}\ttxnova\texon\t{s}\t{e}\t.\t{strand}\t.\t{tattr}")
    dest.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")


def write_universe_gtf(merged: Path, residual_gtf: Path, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w", encoding="utf-8") as out:
        if merged.is_file() and merged.stat().st_size:
            last = ""
            with merged.open(encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    if not is_model_gtf_line(line):
                        continue
                    if 'gene_id "RSDL.' in line or "gene_id 'RSDL." in line:
                        raise TxNovaError(
                            f"{merged} already contains RSDL. gene_id; refuse to concatenate"
                        )
                    out.write(line)
                    last = line
            if last and not last.endswith("\n"):
                out.write("\n")
        if residual_gtf.is_file() and residual_gtf.stat().st_size:
            with residual_gtf.open(encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    out.write(line)


def write_residuals(
    leak_tsv: Path,
    dest: Path,
    *,
    gene_gtfs: list[Path] | None = None,
    samples_json: str | None = None,
    min_mapq: int = 10,
    skip_duplicate: str = "auto",
) -> pd.DataFrame:
    if not leak_tsv.is_file() or leak_tsv.stat().st_size == 0:
        return _write_empty(dest)
    leak = pd.read_csv(leak_tsv, sep="\t")
    if leak.empty:
        return _write_empty(dest)
    genes = merge_gene_bodies(gene_gtfs)
    out, n_deg = cluster_leak(leak, genes=genes or None)
    dest.parent.mkdir(parents=True, exist_ok=True)
    if out.empty:
        _write_empty(dest)
        _write_residual_stats(dest, 0, n_deg)
        return out
    if samples_json:
        out, n_cov = _extend_with_coverage(
            out,
            samples_json,
            genes=genes or None,
            min_mapq=min_mapq,
            skip_duplicate=skip_duplicate,
        )
        n_deg += n_cov
    if out.empty:
        _write_empty(dest)
        _write_residual_stats(dest, 0, n_deg)
        return out
    out.to_csv(dest, sep="\t", index=False)
    write_residual_gtf(out, dest.with_name(RESIDUAL_GTF))
    n_multi = int((out["n_junctions"] >= 2).sum())
    _write_residual_stats(dest, len(out), n_deg)
    log.info(
        "residual loci: %s (%s with ≥2 junctions, %s degenerate dropped)",
        len(out),
        n_multi,
        n_deg,
    )
    return out


def _write_residual_stats(dest: Path, n_loci: int, n_degenerate: int) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    (dest.parent / RESIDUAL_STATS).write_text(
        json.dumps({"n_loci": int(n_loci), "n_degenerate": int(n_degenerate)}) + "\n",
        encoding="utf-8",
    )


def _extend_with_coverage(
    out: pd.DataFrame,
    samples_json: str,
    *,
    genes: GeneBodies | None,
    min_mapq: int,
    skip_duplicate: str,
) -> tuple[pd.DataFrame, int]:
    from txnova import _core

    with tempfile.TemporaryDirectory(prefix="txnova_term_") as td:
        tmp = Path(td)
        loci = tmp / "residual.introns.tsv"
        extents_path = tmp / "residual.terminals.tsv"
        rows: list[str] = ["residual_id\tchrom\tstrand\tintron_start\tintron_end"]
        for rec in out.to_dict(orient="records"):
            introns = parse_exons(rec.get("intron_structure"))
            if not introns:
                continue
            rows.append(
                f"{rec['residual_id']}\t{rec['chrom']}\t{rec['strand']}\t"
                f"{introns[0][0]}\t{introns[-1][1]}"
            )
        loci.write_text("\n".join(rows) + "\n", encoding="utf-8")
        cfg = {
            "min_mapq": min_mapq,
            "skip_duplicate": skip_duplicate,
            "max_terminal_nt": MAX_TERMINAL_NT,
            "min_depth": MIN_TERMINAL_DEPTH,
            "max_gap_nt": MAX_TERMINAL_GAP_NT,
        }
        _core.extend_terminals(str(loci), str(extents_path), samples_json, json.dumps(cfg))
        if not extents_path.is_file() or extents_path.stat().st_size == 0:
            return out, 0
        extents = pd.read_csv(extents_path, sep="\t")
        if extents.empty:
            return out, 0
        return apply_terminal_extents(out, extents, genes)
