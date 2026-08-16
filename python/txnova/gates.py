from __future__ import annotations

from pathlib import Path
from statistics import median

import pandas as pd

from txnova.config import TxNovaConfig
from txnova.errors import TxNovaError
from txnova.logging import get_logger
from txnova.samples import SampleRow

log = get_logger("txnova.gates")

WORST = {"overlap": 4, "extension": 3, "x": 2, "i": 1, "u": 0}

BASE_CANDIDATE_COLUMNS = [
    "locus_id",
    "locus_coord",
    "chrom",
    "start",
    "end",
    "strand",
    "n_exons",
    "length_nt",
    "exon_structure",
    "class",
    "nearest_gene_id",
    "nearest_gene_name",
    "nearest_distance_bp",
    "nearest_strand",
    "nearest_any_gene_id",
    "nearest_any_gene_name",
    "nearest_any_distance_bp",
    "nearest_any_strand",
    "named_gene_name",
    "named_gene_id",
    "named_gene_type",
    "named_overlap",
    "canonical_splice_fraction",
    "coverage_discontinuity",
    "coverage_valley_mean",
    "coverage_gap_mean",
    "control_max_tpm",
    "treat_median_tpm",
    "treat_n_detected",
    "junction_support",
    "junction_support_min",
    "bridge_read_count",
    "has_bridge",
    "gates_passed",
    "representative_transcript_id",
    "baseMean",
    "log2FoldChange",
    "pvalue",
    "padj",
    "de_status",
    "de_pass",
    "longest_orf_aa",
    "orf_complete",
    "coding_score",
    "fickett_score",
    "coding_label",
]


def per_sample_gate_columns(rows: list[SampleRow]) -> list[str]:
    cols: list[str] = []
    for r in rows:
        cols.extend(
            [
                f"{r.sample_id}_junction_support",
                f"{r.sample_id}_bridge_read_count",
                f"{r.sample_id}_tpm",
                f"{r.sample_id}_count",
            ]
        )
    return cols


def candidate_schema(rows: list[SampleRow]) -> list[str]:
    return BASE_CANDIDATE_COLUMNS + per_sample_gate_columns(rows)


UNNAMED_FILENAME = "candidates.unnamed.tsv"


def write_candidates(cand: pd.DataFrame, path: Path, rows: list[SampleRow]) -> None:
    cols = candidate_schema(rows)
    out = cand.copy()
    for c in cols:
        if c not in out.columns:
            out[c] = (
                "NA"
                if c in DE_COLUMNS or c in CODING_COLUMNS
                else ("none" if c == "named_overlap" else pd.Series(dtype=object))
            )
    extras = [c for c in out.columns if c not in cols]
    path.parent.mkdir(parents=True, exist_ok=True)
    out[cols + extras].to_csv(path, sep="\t", index=False, na_rep="NA")


DE_COLUMNS = [
    "baseMean",
    "log2FoldChange",
    "pvalue",
    "padj",
    "de_status",
    "de_pass",
]
CODING_COLUMNS = [
    "longest_orf_aa",
    "orf_complete",
    "coding_score",
    "fickett_score",
    "coding_label",
]

NAMING_COLUMNS = [
    "named_gene_name",
    "named_gene_id",
    "named_gene_type",
    "named_overlap",
]


def drop_join_columns(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    have = [c for c in cols if c in df.columns]
    return df.drop(columns=have) if have else df


STRUCTURE_FEATURE_COLUMNS = (
    "locus_id",
    "sample_id",
    "transcript_id",
    "chrom",
    "start",
    "end",
    "strand",
    "n_exons",
    "length_nt",
    "donors",
    "acceptors",
    "canonical_splice_fraction",
    "n_introns",
    "has_nearest",
    "nearest_gene_id",
    "nearest_gene_name",
    "nearest_distance_bp",
    "nearest_strand",
    "nearest_any_gene_id",
    "nearest_any_gene_name",
    "nearest_any_distance_bp",
    "nearest_any_strand",
    "locus_mean_depth",
    "valley_mean",
    "valley_possible",
    "gap_mean_depth",
    "n_dup_flag_seen",
    "junction_support",
    "bridge_read_count",
    "structure_error",
)


def _csv_fields(value: object) -> list[str]:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return []
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return []
    return [p.strip() for p in text.split(",") if p.strip()]


def is_structure_error(value: object) -> bool:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return False
    s = str(value).strip()
    return bool(s) and s.lower() != "nan"


def validate_structure_features(feat: pd.DataFrame, path: Path) -> None:
    missing = [c for c in STRUCTURE_FEATURE_COLUMNS if c not in feat.columns]
    if missing:
        raise TxNovaError(
            f"{path} missing {', '.join(missing)}; delete stamps/structure.json and rerun"
        )
    for i, rec in enumerate(feat.to_dict(orient="records"), start=2):
        lid = rec.get("locus_id", "")
        sid = rec.get("sample_id", "")
        try:
            n_introns = int(rec["n_introns"])
        except (TypeError, ValueError, KeyError) as e:
            raise TxNovaError(f"{path} line {i} locus {lid}: bad n_introns") from e
        if n_introns < 0:
            raise TxNovaError(f"{path} line {i} locus {lid}: n_introns={n_introns}")
        donors = _csv_fields(rec.get("donors", ""))
        accs = _csv_fields(rec.get("acceptors", ""))
        junc = _csv_fields(rec.get("junction_support", ""))
        if len(donors) != n_introns or len(accs) != n_introns:
            raise TxNovaError(
                f"{path} line {i} locus {lid}: len(donors)={len(donors)} "
                f"len(acceptors)={len(accs)} n_introns={n_introns}"
            )
        if junc:
            if not all(x.isdigit() for x in junc):
                raise TxNovaError(
                    f"{path} line {i} locus {lid} sample {sid}: bad junction_support {junc!r}"
                )
            if len(junc) != n_introns:
                raise TxNovaError(
                    f"{path} line {i} locus {lid} sample {sid}: "
                    f"len(junction_support)={len(junc)} n_introns={n_introns}"
                )


def _median(xs: list[float]) -> float:
    return float(median(xs)) if xs else 0.0


def _n_introns(head: pd.Series | None, n_exons: int) -> int:
    if head is not None:
        raw = head.get("n_introns", "")
        if raw != "" and pd.notna(raw):
            return int(float(raw))
    return max(0, n_exons - 1)


def _splice_fraction(head: pd.Series | None) -> float | None:
    if head is None:
        return None
    raw = head.get("canonical_splice_fraction", "")
    if raw == "" or pd.isna(raw):
        return None
    return float(raw)


def pick_representative_transcript(class_tsv: Path, out_tsv: Path) -> pd.DataFrame:
    """locus.class = u iff every isoform is u. Representative only for those loci."""
    df = pd.read_csv(class_tsv, sep="\t")
    if df.empty:
        out_tsv.parent.mkdir(parents=True, exist_ok=True)
        out_tsv.write_text("locus_id\ttranscript_id\n", encoding="utf-8")
        return pd.DataFrame(columns=["locus_id", "transcript_id"])

    def worst(s: pd.Series) -> str:
        return max(s, key=lambda c: WORST.get(str(c), 0))

    loc = (
        df.groupby("gene_id", sort=False)
        .agg(
            locus_class=("class", worst),
            n_isoforms=("transcript_id", "size"),
            all_u=("class", lambda s: bool((s == "u").all())),
        )
        .reset_index()
        .rename(columns={"gene_id": "locus_id"})
    )
    u_ids = set(loc.loc[loc["all_u"], "locus_id"])
    cand = df[df["gene_id"].isin(u_ids)].copy()
    if cand.empty:
        out_tsv.parent.mkdir(parents=True, exist_ok=True)
        out_tsv.write_text("locus_id\ttranscript_id\n", encoding="utf-8")
        loc.to_csv(out_tsv.with_name("locus.class.tsv"), sep="\t", index=False)
        return pd.DataFrame(columns=["locus_id", "transcript_id"])

    cand["length_nt"] = pd.to_numeric(cand["length_nt"], errors="coerce").fillna(0)
    cand = cand.sort_values(
        ["gene_id", "length_nt", "transcript_id"], ascending=[True, False, True]
    )
    reps = cand.groupby("gene_id", sort=False).first().reset_index()
    out = reps.rename(columns={"gene_id": "locus_id"})[["locus_id", "transcript_id"]]
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_tsv, sep="\t", index=False)
    loc.to_csv(out_tsv.with_name("locus.class.tsv"), sep="\t", index=False)
    return out


def _locus_record(
    *,
    t: pd.Series,
    lid: str,
    tid: str,
    n_exons: int,
    length_nt: int,
    head: pd.Series | None,
    frac: float,
    disc_ok: bool,
    valley_med: object,
    gap_med: object,
    control_max: float,
    treat_med: float,
    n_det: int,
    junc_sum: list[int] | None,
    bridge_total: int,
    has_bridge: bool,
    passed: list[str],
    sample_junc: dict[str, str],
    sample_bridge: dict[str, int],
    tpm_row: pd.DataFrame,
    counts: pd.DataFrame,
    rows: list[SampleRow],
) -> dict:
    chrom, start, end, strand = t["chrom"], int(t["start"]), int(t["end"]), t["strand"]
    locus_coord = f"{chrom}:{start}-{end}:{strand}"
    exons = "" if "exon_structure" not in t.index else str(t["exon_structure"])
    if exons == "nan":
        exons = ""
    rec_out: dict = {
        "locus_id": lid,
        "locus_coord": locus_coord,
        "chrom": chrom,
        "start": start,
        "end": end,
        "strand": strand,
        "n_exons": n_exons,
        "length_nt": length_nt,
        "exon_structure": exons,
        "class": "u",
        "nearest_gene_id": "" if head is None else head.get("nearest_gene_id", ""),
        "nearest_gene_name": "" if head is None else head.get("nearest_gene_name", ""),
        "nearest_distance_bp": "" if head is None else head.get("nearest_distance_bp", ""),
        "nearest_strand": "" if head is None else head.get("nearest_strand", ""),
        "nearest_any_gene_id": "" if head is None else head.get("nearest_any_gene_id", ""),
        "nearest_any_gene_name": "" if head is None else head.get("nearest_any_gene_name", ""),
        "nearest_any_distance_bp": "" if head is None else head.get("nearest_any_distance_bp", ""),
        "nearest_any_strand": "" if head is None else head.get("nearest_any_strand", ""),
        "canonical_splice_fraction": frac,
        "coverage_discontinuity": disc_ok,
        "coverage_valley_mean": valley_med,
        "coverage_gap_mean": gap_med,
        "control_max_tpm": control_max,
        "treat_median_tpm": treat_med,
        "treat_n_detected": n_det,
        "junction_support": ",".join(str(x) for x in junc_sum) if junc_sum else "",
        "junction_support_min": min(junc_sum) if junc_sum else "",
        "bridge_read_count": bridge_total,
        "has_bridge": has_bridge,
        "gates_passed": ",".join(passed),
        "representative_transcript_id": tid,
    }
    for sid, js in sample_junc.items():
        rec_out[f"{sid}_junction_support"] = js
        rec_out[f"{sid}_bridge_read_count"] = sample_bridge.get(sid, 0)
    if not tpm_row.empty:
        for sid in [r.sample_id for r in rows]:
            rec_out[f"{sid}_tpm"] = float(tpm_row.iloc[0][sid]) if sid in tpm_row.columns else 0.0
    if not counts.empty:
        crow = counts[counts["locus_id"] == lid]
        for sid in [r.sample_id for r in rows]:
            rec_out[f"{sid}_count"] = (
                float(crow.iloc[0][sid]) if not crow.empty and sid in crow.columns else 0.0
            )
    return rec_out


def apply_gates(
    cfg: TxNovaConfig,
    rows: list[SampleRow],
    *,
    class_tsv: Path,
    reps_tsv: Path,
    structure_tsv: Path,
    locus_tpm: Path,
    locus_counts: Path,
    out_dir: Path,
    filename: str = "candidates.tsv",
) -> pd.DataFrame:
    if not reps_tsv.is_file():
        raise TxNovaError(f"missing {reps_tsv}")
    reps = pd.read_csv(reps_tsv, sep="\t")
    cls = pd.read_csv(class_tsv, sep="\t")
    feat = pd.read_csv(structure_tsv, sep="\t") if structure_tsv.is_file() else pd.DataFrame()
    if not feat.empty:
        validate_structure_features(feat, structure_tsv)
    tpm = pd.read_csv(locus_tpm, sep="\t") if locus_tpm.is_file() else pd.DataFrame()
    counts = pd.read_csv(locus_counts, sep="\t") if locus_counts.is_file() else pd.DataFrame()
    if not reps.empty:
        if tpm.empty or "locus_id" not in tpm.columns:
            raise TxNovaError(f"missing or empty full-universe TPM table {locus_tpm}")
        have = set(tpm["locus_id"].astype(str))
        need = {str(x) for x in reps["locus_id"]}
        missing = sorted(need - have)
        if missing:
            raise TxNovaError(
                f"{locus_tpm} missing locus_id {missing[0]}"
                + (f" and {len(missing) - 1} others" if len(missing) > 1 else "")
                + "; full-universe TPM is required"
            )

    controls = [r.sample_id for r in rows if r.group == "control"]
    treats = [r.sample_id for r in rows if r.group == "treat"]
    f = cfg.filters
    out_rows: list[dict] = []
    unnamed_rows: list[dict] = []
    uns = all(r.strandedness == "unstranded" for r in rows)

    for rec in reps.to_dict(orient="records"):
        lid = rec["locus_id"]
        tid = rec["transcript_id"]
        iso = cls[cls["transcript_id"] == tid]
        if iso.empty:
            continue
        t = iso.iloc[0]
        sf = feat[feat["locus_id"] == lid] if not feat.empty else pd.DataFrame()
        head = sf.iloc[0] if not sf.empty else None
        if head is not None and is_structure_error(head.get("structure_error", "")):
            log.warning("locus %s skipped: %s", lid, head.get("structure_error", ""))
            continue

        passed = ["class"]
        n_exons = int(t["n_exons"])
        if n_exons < f.min_exons:
            continue
        passed.append("min_exons")
        length_nt = int(t["length_nt"])
        if length_nt < f.transcript_min_nt:
            continue
        passed.append("transcript_min_nt")

        n_introns = _n_introns(head, n_exons)
        frac = _splice_fraction(head)
        if n_introns > 0:
            scored = 0.0 if frac is None else frac
            if f.require_canonical_splice and scored < (1.0 - f.max_noncanonical_junction_fraction):
                continue
            passed.append("canonical_splice")
        elif frac is None:
            frac = ""

        has_nearest = False
        dist = None
        if head is not None:
            has_nearest = str(head.get("has_nearest", "")).lower() == "true"
            col = "nearest_any_distance_bp" if uns else "nearest_distance_bp"
            raw = head.get(col, "")
            if raw != "" and pd.notna(raw):
                dist = float(raw)
        if has_nearest and dist is None:
            raise TxNovaError(f"locus {lid}: has_nearest=true but nearest distance is empty")
        if has_nearest and dist is not None and dist < f.min_nearest_same_strand_bp:
            continue
        passed.append("distance")

        tpm_row = tpm[tpm["locus_id"] == lid] if not tpm.empty else pd.DataFrame()
        treat_tpms_by = {
            sid: (
                float(tpm_row.iloc[0][sid]) if not tpm_row.empty and sid in tpm_row.columns else 0.0
            )
            for sid in treats
        }
        detecting = [sid for sid, tv in treat_tpms_by.items() if tv >= f.treat_detect_tpm]

        # discontinuity
        disc_ok = True
        valley_med = ""
        gap_med = ""
        if not has_nearest:
            passed.append("discontinuity")
        elif f.require_coverage_discontinuity:
            if sf.empty or not detecting:
                disc_ok = False
            else:
                ok_n = 0
                detecting_valleys: list[float] = []
                detecting_gaps: list[float] = []
                for sid in detecting:
                    sub = sf[sf["sample_id"] == sid]
                    if sub.empty:
                        continue
                    row = sub.iloc[0]
                    vp = str(row.get("valley_possible", "")).lower()
                    if vp != "true":
                        continue
                    try:
                        vm = float(row["valley_mean"])
                        lm = float(row["locus_mean_depth"])
                        gm = float(row["gap_mean_depth"])
                    except (TypeError, ValueError):
                        continue
                    detecting_valleys.append(vm)
                    detecting_gaps.append(gm)
                    thresh = min(f.discontinuity_valley_max_mean, f.discontinuity_ratio * lm)
                    # valley: a local dip exists. gap_mean: the intervening
                    # sequence is actually empty, not a transcribed 3′ tail
                    # with one empty 200 bp window.
                    if vm <= thresh and gm <= thresh:
                        ok_n += 1
                need = (
                    len(detecting)
                    if f.discontinuity_min_treat_samples is None
                    else f.discontinuity_min_treat_samples
                )
                if ok_n < need:
                    disc_ok = False
                if detecting_valleys:
                    valley_med = _median(detecting_valleys)
                if detecting_gaps:
                    gap_med = _median(detecting_gaps)
            if not disc_ok:
                continue
            passed.append("discontinuity")
        else:
            passed.append("discontinuity")

        bridge_total = 0
        junc_sum: list[int] | None = None
        sample_junc: dict[str, str] = {}
        sample_bridge: dict[str, int] = {}
        if not sf.empty:
            for sid in [r.sample_id for r in rows]:
                sub = sf[sf["sample_id"] == sid]
                if sub.empty:
                    sample_junc[sid] = ""
                    sample_bridge[sid] = 0
                    continue
                row = sub.iloc[0]
                js = (
                    ""
                    if pd.isna(row.get("junction_support", ""))
                    else str(row.get("junction_support", ""))
                )
                if js == "nan":
                    js = ""
                sample_junc[sid] = js
                try:
                    br = int(float(row.get("bridge_read_count", 0) or 0))
                except (TypeError, ValueError):
                    br = 0
                sample_bridge[sid] = br
            for sid in detecting:
                bridge_total += sample_bridge.get(sid, 0)
                fields = _csv_fields(sample_junc.get(sid, ""))
                if not fields:
                    continue
                parts = [int(x) for x in fields]
                if junc_sum is None:
                    junc_sum = parts
                else:
                    junc_sum = [a + b for a, b in zip(junc_sum, parts, strict=True)]
        has_bridge = bool(has_nearest and detecting and bridge_total >= f.bridge_min_reads)
        if f.reject_bridging_junction:
            if has_bridge:
                continue
            passed.append("no_bridge")

        # detection
        ctrl_tpms = [
            float(tpm_row.iloc[0][sid]) if not tpm_row.empty and sid in tpm_row.columns else 0.0
            for sid in controls
        ]
        treat_tpms = [treat_tpms_by[sid] for sid in treats]
        control_max = max(ctrl_tpms) if ctrl_tpms else 0.0
        treat_med = _median(treat_tpms)
        n_det = len(detecting)
        rec_out = _locus_record(
            t=t,
            lid=lid,
            tid=tid,
            n_exons=n_exons,
            length_nt=length_nt,
            head=head,
            frac=frac,
            disc_ok=disc_ok,
            valley_med=valley_med,
            gap_med=gap_med,
            control_max=control_max,
            treat_med=treat_med,
            n_det=n_det,
            junc_sum=junc_sum,
            bridge_total=bridge_total,
            has_bridge=has_bridge,
            passed=passed,
            sample_junc=sample_junc,
            sample_bridge=sample_bridge,
            tpm_row=tpm_row,
            counts=counts,
            rows=rows,
        )
        if control_max >= f.control_max_tpm:
            unnamed_rows.append(rec_out)
            continue
        if n_det < f.treat_min_detected_replicates:
            continue
        if treat_med < f.treat_median_tpm:
            continue
        rec_out["gates_passed"] = ",".join(passed + ["detect"])
        out_rows.append(rec_out)

    out_dir.mkdir(parents=True, exist_ok=True)
    cand = pd.DataFrame(out_rows)
    unnamed = pd.DataFrame(unnamed_rows)
    write_candidates(cand, out_dir / filename, rows)
    write_candidates(unnamed, out_dir / UNNAMED_FILENAME, rows)
    log.info("candidates: %s treat-specific, %s unnamed (also in control)", len(cand), len(unnamed))
    return cand
