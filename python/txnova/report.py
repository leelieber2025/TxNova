from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd
from jinja2 import Environment, FileSystemLoader, select_autoescape

from txnova.config import TxNovaConfig
from txnova.gates import is_structure_error
from txnova.samples import SampleRow
from txnova.staging import StagingDir

_TPL = Path(__file__).resolve().parent / "templates"

CANDIDATE_COLS = [
    "locus_id",
    "locus_coord",
    "n_exons",
    "nearest_gene_name",
    "named_gene_name",
    "named_overlap",
    "nearest_distance_bp",
    "longest_orf_aa",
    "coverage_valley_mean",
    "coverage_gap_mean",
    "assembled_in",
    "junction_support",
    "bridge_read_count",
    "control_max_tpm",
    "treat_median_tpm",
    "baseMean",
    "pvalue",
    "padj",
    "de_status",
    "coding_label",
]


def _fmt_bytes(n: Any) -> str:
    if not isinstance(n, (int, float)) or n <= 0:
        return "unknown"
    return f"{n / (1024**3):.1f} GiB"


def _threads_line(preflight: dict[str, Any]) -> str:
    t = preflight.get("threads") or {}
    return (
        f"threads: BAM ceiling {t.get('bam_workers')} "
        f"(live follows MemAvailable; cpus {t.get('cpus')}, usable {t.get('usable')}, "
        f"MemAvailable {_fmt_bytes(t.get('available_bytes'))}; "
        f"requested={t.get('requested')}; auto={t.get('auto')})"
    )


def _cell(v: Any) -> str:
    if v is None or (isinstance(v, float) and pd.isna(v)) or (isinstance(v, str) and v == "nan"):
        return ""
    return str(v)


def _funnel(out: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    cls = out / "classify" / "transcripts.class.tsv"
    loc = out / "classify" / "locus.class.tsv"
    cand = out / "candidates" / "candidates.tsv"
    if cls.is_file():
        df = pd.read_csv(cls, sep="\t")
        rows.append({"label": "merged transcripts", "n": int(len(df))})
        if "class" in df.columns:
            rows.append({"label": "class u transcripts", "n": int((df["class"] == "u").sum())})
    if loc.is_file():
        df = pd.read_csv(loc, sep="\t")
        rows.append({"label": "merged loci", "n": int(len(df))})
        if "all_u" in df.columns:
            rows.append(
                {
                    "label": "all-u loci",
                    "n": int(df["all_u"].astype(str).str.lower().eq("true").sum()),
                }
            )
        elif "locus_class" in df.columns:
            rows.append({"label": "all-u loci", "n": int((df["locus_class"] == "u").sum())})
    struct = out / "structure" / "structure.features.tsv"
    if struct.is_file():
        sdf = pd.read_csv(struct, sep="\t")
        if not sdf.empty and "structure_error" in sdf.columns and "locus_id" in sdf.columns:
            n_err = int(
                sdf.loc[sdf["structure_error"].map(is_structure_error), "locus_id"].nunique()
            )
            rows.append({"label": "structure_error loci", "n": n_err})
    gates = out / "candidates" / "candidates.gates.tsv"
    unnamed = out / "candidates" / "candidates.unnamed.tsv"
    de_view = out / "candidates" / "candidates.de.tsv"
    if unnamed.is_file():
        rows.append(
            {"label": "unnamed (structure-pass, also in control)", "n": _tsv_nrows(unnamed)}
        )
    if gates.is_file():
        rows.append({"label": "after gates", "n": _tsv_nrows(gates)})
    if de_view.is_file():
        rows.append({"label": "after DE", "n": _tsv_nrows(de_view)})
    if cand.is_file():
        rows.append({"label": "final candidates", "n": _tsv_nrows(cand)})
    orphan = out / "candidates" / "orphan.tsv"
    if orphan.is_file():
        rows.append({"label": "orphan (no gene name)", "n": _tsv_nrows(orphan)})
    leak = out / "candidates" / "leak.tsv"
    if leak.is_file():
        ldf = pd.read_csv(leak, sep="\t")
        if ldf.empty or "status" not in ldf.columns:
            rows.append({"label": "leak (treat-specific splice, not a candidate)", "n": 0})
        else:
            n_un = int((ldf["status"] == "unassembled").sum())
            n_au = int((ldf["status"] == "assembled_u").sum())
            rows.append({"label": "leak unassembled (BAM splice, not in annotation)", "n": n_un})
            if n_au:
                rows.append(
                    {
                        "label": "leak assembled_u (unused; merged.gtf is the annotation)",
                        "n": n_au,
                    }
                )
    residual = out / "candidates" / "residual.tsv"
    if residual.is_file():
        rdf = _read_tsv(residual)
        if rdf.empty:
            rows.append({"label": "residual splice loci", "n": 0})
        else:
            n_all = int(len(rdf))
            n_multi = int((rdf["n_junctions"] >= 2).sum()) if "n_junctions" in rdf.columns else 0
            n_geneish = 0
            if {"n_exons", "treat_n_detected"}.issubset(rdf.columns):
                n_geneish = int(((rdf["n_exons"] >= 3) & (rdf["treat_n_detected"] >= 3)).sum())
            rows.append({"label": "residual splice loci (unassembled intergenic)", "n": n_all})
            rows.append({"label": "residual ≥2 junctions", "n": n_multi})
            rows.append({"label": "residual ≥3 exons, treat_n≥3", "n": n_geneish})
    gr = out / "candidates" / "gene_rank.tsv"
    if gr.is_file():
        gdf = _read_tsv(gr)
        rows.append(
            {"label": "gene-like ranked (structure-pass)", "n": 0 if gdf.empty else int(len(gdf))}
        )
    return rows


def _tsv_nrows(path: Path) -> int:
    return max(0, sum(1 for _ in path.open()) - 1)


LEAK_COLS = [
    "chrom",
    "start",
    "end",
    "strand",
    "status",
    "merged_locus",
    "merged_class",
    "overlaps_gene",
    "nearest_gene_name",
    "nearest_distance_bp",
    "control_max",
    "treat_sum",
    "treat_n_detected",
    "in_gates",
    "in_unnamed",
]


def _leak_records(path: Path) -> tuple[list[str], list[dict[str, str]], int]:
    if not path.is_file():
        return [], [], 0
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return [], [], 0
    n = int(len(df))
    if "treat_sum" in df.columns:
        df = df.sort_values("treat_sum", ascending=False, na_position="last")
    df = df.head(40)
    cols = [c for c in LEAK_COLS if c in df.columns]
    recs = [{c: _cell(row[c]) for c in cols} for row in df.to_dict(orient="records")]
    return cols, recs, n


RESIDUAL_COLS = [
    "residual_id",
    "locus_coord",
    "n_junctions",
    "n_exons",
    "length_nt",
    "nearest_gene_name",
    "nearest_distance_bp",
    "treat_sum",
    "treat_n_detected",
    "n_shared_site",
]


def _read_tsv(path: Path) -> pd.DataFrame:
    if not path.is_file() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t")
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def _residual_records(path: Path) -> tuple[list[str], list[dict[str, str]], int]:
    df = _read_tsv(path)
    if df.empty:
        return [], [], 0
    n = int(len(df))
    sort_cols = [c for c in ("n_junctions", "treat_sum", "treat_n_detected") if c in df.columns]
    if sort_cols:
        df = df.sort_values(sort_cols, ascending=False, na_position="last")
    df = df.head(40)
    cols = [c for c in RESIDUAL_COLS if c in df.columns]
    recs = [{c: _cell(row[c]) for c in cols} for row in df.to_dict(orient="records")]
    return cols, recs, n


GENE_RANK_COLS = [
    "gene_rank",
    "locus_id",
    "locus_coord",
    "n_exons",
    "length_nt",
    "nearest_gene_name",
    "nearest_distance_bp",
    "treat_median_tpm",
    "control_max_tpm",
    "longest_orf_aa",
    "coding_score",
    "coding_label",
    "gene_score",
    "gene_flags",
    "phylop_mean",
    "pfam_name",
    "pfam_evalue",
    "pfam_desc",
]


def _gene_rank_records(
    path: Path, *, limit: int = 30
) -> tuple[list[str], list[dict[str, str]], int]:
    df = _read_tsv(path)
    if df.empty:
        return [], [], 0
    n = int(len(df))
    if "gene_rank" in df.columns:
        df = df.sort_values("gene_rank", ascending=True)
    df = df.head(limit)
    cols = [c for c in GENE_RANK_COLS if c in df.columns]
    recs = [{c: _cell(row[c]) for c in cols} for row in df.to_dict(orient="records")]
    return cols, recs, n


ORPHAN_COLS = [
    "locus_id",
    "chrom",
    "longest_orf_aa",
    "coding_label",
    "phylop_mean",
    "phylop_frac_pos",
    "phastcons_mean",
    "pfam_name",
    "pfam_evalue",
    "pfam_desc",
]


def _orphan_records(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file():
        return [], []
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return [], []
    cols = [c for c in ORPHAN_COLS if c in df.columns]
    recs = [{c: _cell(row[c]) for c in cols} for row in df.to_dict(orient="records")]
    return cols, recs


def _candidate_records(
    path: Path,
    *,
    limit: int | None = None,
    sort_col: str | None = None,
) -> tuple[list[str], list[dict[str, str]], int]:
    if not path.is_file():
        return [], [], 0
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return [], [], 0
    n = int(len(df))
    if sort_col and sort_col in df.columns:
        df = df.sort_values(sort_col, ascending=False, na_position="last")
    if limit is not None:
        df = df.head(limit)
    cols = [c for c in CANDIDATE_COLS if c in df.columns]
    recs = [{c: _cell(row[c]) for c in cols} for row in df.to_dict(orient="records")]
    return cols, recs, n


def build_context(
    cfg: TxNovaConfig,
    rows: list[SampleRow],
    preflight: dict[str, Any],
    n_candidates: int = 0,
) -> dict[str, Any]:
    cand_path = cfg.output_dir / "candidates" / "candidates.tsv"
    unnamed_path = cfg.output_dir / "candidates" / "candidates.unnamed.tsv"
    cols, recs, n_from_file = _candidate_records(cand_path)
    if cand_path.is_file():
        n_candidates = n_from_file
    u_cols, u_recs, n_unnamed = _candidate_records(unnamed_path)
    orphan_path = cfg.output_dir / "candidates" / "orphan.tsv"
    o_cols, o_recs = _orphan_records(orphan_path)
    leak_path = cfg.output_dir / "candidates" / "leak.tsv"
    l_cols, l_recs, n_leaks = _leak_records(leak_path)
    residual_path = cfg.output_dir / "candidates" / "residual.tsv"
    r_cols, r_recs, n_residuals = _residual_records(residual_path)
    rank_path = cfg.output_dir / "candidates" / "gene_rank.tsv"
    g_cols, g_recs, n_ranked = _gene_rank_records(rank_path, limit=30)
    return {
        "cfg": cfg,
        "rows": rows,
        "preflight": preflight,
        "n_candidates": n_candidates,
        "warnings": preflight.get("warnings") or [],
        "errors": preflight.get("errors") or [],
        "funnel": _funnel(cfg.output_dir),
        "candidate_cols": cols,
        "candidates": recs,
        "unnamed_cols": u_cols,
        "unnamed": u_recs,
        "n_unnamed": n_unnamed,
        "has_structures": (cfg.output_dir / "fold" / "index.html").is_file(),
        "orphan_cols": o_cols,
        "orphans": o_recs,
        "n_orphans": len(o_recs),
        "leak_cols": l_cols,
        "leaks": l_recs,
        "n_leaks": n_leaks,
        "residual_cols": r_cols,
        "residuals": r_recs,
        "n_residuals": n_residuals,
        "gene_rank_cols": g_cols,
        "gene_ranks": g_recs,
        "n_ranked": n_ranked,
    }


def render_report(
    cfg: TxNovaConfig,
    rows: list[SampleRow],
    preflight: dict[str, Any],
    n_candidates: int = 0,
) -> str:
    ctx = build_context(cfg, rows, preflight, n_candidates=n_candidates)
    lines = [
        "# TxNova report",
        "",
        f"**{ctx['n_candidates']}** locus/loci in the final table.",
        "",
        "## Run",
        "",
        f"- species: `{cfg.species}`",
        f"- assembly: `{cfg.genome.assembly}`",
        f"- annotation: `{cfg.genome.annotation_source} {cfg.genome.annotation_version}`",
        f"- aligner family: `{preflight.get('aligner_family', '')}`",
        f"- library layout: `{preflight.get('library_layout', '')}`",
        f"- strandedness: `{preflight.get('strandedness', '')}`",
        f"- control / treat: {preflight.get('n_control', 0)} / {preflight.get('n_treat', 0)}",
        f"- {_threads_line(preflight)}",
        f"- DE: {'DESeq2-style (pydeseq2); keeps wald (padj+LFC) and low_count (independent-filtering NA); evaluated non-hits are dropped' if cfg.de.enabled else 'not run'}",
        "- coding: " + ("hexamer LLR + Fickett" if cfg.coding.enabled else "not run"),
        "- orphan: "
        + (
            "named_overlap=none: UCSC phyloP/phastCons + HMMER hmmscan/Pfam"
            if cfg.coding.enabled and cfg.coding.orphan
            else "not run"
        ),
        "",
        "## Samples",
        "",
        "| sample_id | group | replicate | strandedness | bam |",
        "| --- | --- | --- | --- | --- |",
    ]
    for r in rows:
        lines.append(
            f"| {r.sample_id} | {r.group} | {r.replicate} | {r.strandedness} | `{r.bam}` |"
        )
    lines += [
        "",
        "## Preflight",
        "",
        f"- ok: `{preflight.get('ok')}`",
        f"- SQ sha256: `{preflight.get('sq_sha256', '')}`",
        f"- GTF seqnames: {preflight.get('gtf_seqnames', '')}",
        "",
    ]
    if ctx["warnings"]:
        lines.append("### Warnings")
        lines.append("")
        for w in ctx["warnings"]:
            lines.append(f"- {w}")
        lines.append("")
    if ctx["errors"]:
        lines.append("### Errors")
        lines.append("")
        for e in ctx["errors"]:
            lines.append(f"- {e}")
        lines.append("")
    if ctx["funnel"]:
        lines += ["## Funnel", "", "| step | n |", "| --- | --- |"]
        for row in ctx["funnel"]:
            lines.append(f"| {row['label']} | {row['n']} |")
        lines += [
            "",
            "Annotation plus residual splice loci form the quantify universe. Residual loci cluster treat-specific "
            "BAM splices missing from the annotation. Class `u` is a negative "
            "definition. Counts and junctions are recomputed from the BAM.",
            "",
        ]
    if ctx.get("has_structures"):
        lines += [
            "## Structures",
            "",
            "3D models for ORFs: named loci from AlphaFold DB, unnamed peptides from ESMFold. "
            "Colour is pLDDT (blue = confident). See `fold/index.html`.",
            "",
        ]
        fn = cfg.output_dir / "fold" / "function.tsv"
        if fn.is_file():
            import pandas as pd

            fdf = pd.read_csv(fn, sep="\t")
            if not fdf.empty:
                lines += [
                    "Function from UniProt (if named) and Foldseek fold matches "
                    "(`fold/function.tsv`). A Foldseek hit is a similar 3D fold, "
                    "not a measured activity.",
                    "",
                ]
                show = [
                    c
                    for c in (
                        "locus_id",
                        "gene",
                        "curated_function",
                        "fold_description",
                        "fold_evalue",
                    )
                    if c in fdf.columns
                ]
                lines.append("| " + " | ".join(show) + " |")
                lines.append("| " + " | ".join("---" for _ in show) + " |")
                for rec in fdf.to_dict(orient="records"):
                    cells = []
                    for c in show:
                        v = rec.get(c, "")
                        if isinstance(v, float) and pd.isna(v):
                            v = ""
                        s = str(v).replace("|", "\\|")
                        if len(s) > 160:
                            s = s[:157] + "..."
                        cells.append(s)
                    lines.append("| " + " | ".join(cells) + " |")
                lines.append("")
    lines += [
        "## Leak (BAM splice not in the candidate path)",
        "",
    ]
    if ctx.get("leaks"):
        cols = ctx["leak_cols"]
        shown = len(ctx["leaks"])
        extra = f" Showing {shown} with highest treat support." if shown < ctx["n_leaks"] else ""
        lines.append(
            f"**{ctx['n_leaks']}** treat-recurrent, control-silent splice junctions "
            "missing from the annotation (`unassembled`). Known-gene introns are "
            "omitted. `assembled_u` is unused because `merged.gtf` is the "
            "annotation. This is a recall check, not a gene list."
            f"{extra} See `candidates/leak.tsv`."
        )
        lines.append("")
        lines.append("| " + " | ".join(cols) + " |")
        lines.append("| " + " | ".join("---" for _ in cols) + " |")
        for rec in ctx["leaks"]:
            lines.append("| " + " | ".join(rec[c].replace("|", "\\|") for c in cols) + " |")
        lines.append("")
    else:
        lines += [
            "No treat-specific residual splices. See `candidates/leak.tsv`.",
            "",
        ]
    lines += [
        "## Residual splice loci",
        "",
    ]
    if ctx.get("residuals"):
        cols = ctx["residual_cols"]
        shown = len(ctx["residuals"])
        extra = f" Showing {shown} with most junctions." if shown < ctx["n_residuals"] else ""
        lines.append(
            f"**{ctx['n_residuals']}** intergenic unassembled leak junctions "
            "clustered into locus hypotheses (shared splice site or a 30–20 kb "
            "constitutive exon between adjacent introns). Terminal exons "
            "follow treat coverage from the splice site "
            "(2 kb cap). Residual loci are in `universe.gtf` and go through "
            "quantify, gates, and DE. This table is the locus-hypothesis "
            "list, not the final candidate table."
            f"{extra} See `candidates/residual.tsv`."
        )
        lines.append("")
        lines.append("| " + " | ".join(cols) + " |")
        lines.append("| " + " | ".join("---" for _ in cols) + " |")
        for rec in ctx["residuals"]:
            lines.append("| " + " | ".join(rec[c].replace("|", "\\|") for c in cols) + " |")
        lines.append("")
    else:
        lines += [
            "No residual splice loci. See `candidates/residual.tsv`.",
            "",
        ]
    lines += [
        "## Gene-like rank",
        "",
    ]
    if ctx.get("gene_ranks"):
        cols = ctx["gene_rank_cols"]
        shown = len(ctx["gene_ranks"])
        extra = (
            f" Showing the top {shown} after structure/function on the head of the list."
            if shown < ctx["n_ranked"]
            else ""
        )
        lines.append(
            f"**{ctx['n_ranked']}** structure-pass loci ranked by how gene-like they "
            "are (ORF + hexamer + Fickett; spliced / intergenic; penalties for "
            "unplaced contig, unsupported intron, gag/pol/MLV, intronless copy of "
            "a known protein). phyloP + Pfam on the top 40, then re-sorted."
            f"{extra} See `candidates/gene_rank.tsv`."
        )
        lines.append("")
        lines.append("| " + " | ".join(cols) + " |")
        lines.append("| " + " | ".join("---" for _ in cols) + " |")
        for rec in ctx["gene_ranks"]:
            lines.append("| " + " | ".join(rec[c].replace("|", "\\|") for c in cols) + " |")
        lines.append("")
    else:
        lines += [
            "No gene-like ranking. See `candidates/gene_rank.tsv`.",
            "",
        ]
    if ctx.get("n_unnamed"):
        lines += [
            f"{ctx['n_unnamed']} structure-pass loci are also in control "
            "(`candidates/candidates.unnamed.tsv`). They are included in the rank.",
            "",
        ]
    lines += [
        "## Orphan (no gene name)",
        "",
    ]
    if ctx.get("orphans"):
        cols = ctx["orphan_cols"]
        lines.append(
            f"**{ctx['n_orphans']}** loci with `named_overlap=none`. "
            "phyloP/phastCons are UCSC conservation on exons. "
            "Pfam is HMMER hmmscan. See `candidates/orphan.tsv`."
        )
        lines.append("")
        lines.append("| " + " | ".join(cols) + " |")
        lines.append("| " + " | ".join("---" for _ in cols) + " |")
        for rec in ctx["orphans"]:
            lines.append("| " + " | ".join(rec[c].replace("|", "\\|") for c in cols) + " |")
        lines.append("")
    else:
        lines += [
            "No unnamed-overlap loci were scored. See `candidates/orphan.tsv`.",
            "",
        ]
    lines += [
        "## Candidates",
        "",
    ]
    if ctx["candidates"]:
        cols = ctx["candidate_cols"]
        lines.append("| " + " | ".join(cols) + " |")
        lines.append("| " + " | ".join("---" for _ in cols) + " |")
        for rec in ctx["candidates"]:
            lines.append("| " + " | ".join(rec[c].replace("|", "\\|") for c in cols) + " |")
        lines += [
            "",
            "`assembled_in` is unused with the residual assembler "
            "(per-sample GTFs are stubs). Junction and bridge counts are "
            "BAM CIGAR `N`. `de_status=wald` is padj + LFC. "
            "`de_status=low_count` is independent-filtering NA (Wald was "
            "computed); it is not a significant padj.",
            "",
        ]
    else:
        lines += [
            "Empty final table. Structure + detection may still have rows in "
            "`candidates/candidates.gates.tsv`. Evaluated DE non-hits are "
            "dropped; independent-filtering NA is not. See "
            "`candidates/candidates.tsv` and `quantify/de.tsv`.",
            "",
        ]
    return "\n".join(lines)


def render_html(
    cfg: TxNovaConfig,
    rows: list[SampleRow],
    preflight: dict[str, Any],
    n_candidates: int = 0,
) -> str:
    ctx = build_context(cfg, rows, preflight, n_candidates=n_candidates)
    env = Environment(
        loader=FileSystemLoader(str(_TPL)),
        autoescape=select_autoescape(["html", "j2"]),
    )
    return env.get_template("report.html.j2").render(**ctx)


def write_report(cfg: TxNovaConfig, markdown: str, html: str | None = None) -> Path:
    with StagingDir(cfg.output_dir) as st:
        st.write_text("report/report.md", markdown)
        names = ["report/report.md"]
        if html is not None:
            st.write_text("report/report.html", html)
            names.append("report/report.html")
        st.publish(names)
    return cfg.output_dir / "report" / "report.md"
