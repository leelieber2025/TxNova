from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from txnova.config import TxNovaConfig, resolve_hexamer_table
from txnova.errors import TxNovaError
from txnova.assembly_evidence import attach_assembly_evidence
from txnova.gates import (
    CODING_COLUMNS,
    DE_COLUMNS,
    UNNAMED_FILENAME,
    apply_gates,
    drop_join_columns,
    is_structure_error,
    pick_representative_transcript,
    write_candidates,
)
from txnova.fold import fold_peptides
from txnova.function import annotate_structures
from txnova.naming import label_candidate_file, load_gene_bodies
from txnova.orphan import annotate_orphans
from txnova.leak import LEAK_FILENAME, SHARED_FILENAME, annotate_leak, write_shared
from txnova.residual import (
    CLOSE_SAME_STRAND_NT,
    MAX_EXON_NT,
    MAX_TERMINAL_NT,
    MIN_EXON_NT,
    MIN_INTRON_NT,
    MIN_TERMINAL_NT,
    RESIDUAL_FILENAME,
    RESIDUAL_GTF,
    RESIDUAL_STATS,
    UNIVERSE_GTF,
    write_residuals,
    write_universe_gtf,
)
from txnova.gene_score import GENE_RANK_FILENAME, REPORT_N, top_rank_ids, write_gene_rank
from txnova.logging import get_logger
from txnova.preflight import require_ok, run_preflight, write_preflight_json
from txnova.provenance import build_run_json, write_run_json
from txnova.report import render_html, render_report, write_report
from txnova.samples import SampleRow, can_run_de, load_samples, samples_to_jsonable
from txnova.staging import remove_orphan_staging
from txnova.stamps import (
    bam_fingerprint,
    clear_stamps,
    path_fingerprint,
    sheet_sha,
    stamp_matches,
    stamp_outputs,
    write_stamp,
)
from txnova.assemble import ASSEMBLER_VERSION, Assembler

log = get_logger("txnova.orchestrator")


def _core():
    try:
        from txnova import _core as core
    except ImportError as e:
        raise TxNovaError("txnova._core is not built. pip install -e .") from e
    return core


def _apply_contrast_policy(cfg: TxNovaConfig, rows: list[SampleRow]) -> None:
    if cfg.de.enabled and not can_run_de(rows):
        log.info("DE skipped: need ≥2 control and ≥2 treat")
        cfg.de.enabled = False


def preflight_only(cfg: TxNovaConfig) -> dict[str, Any]:
    rows = load_samples(cfg.samples)
    _apply_contrast_policy(cfg, rows)
    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    remove_orphan_staging(cfg.output_dir)
    report = run_preflight(cfg, rows)
    write_preflight_json(cfg, report)
    require_ok(report)
    return report


def _dup_map(preflight: dict[str, Any]) -> dict[str, bool]:
    raw = preflight.get("dup_flag_seen") or {}
    return {k: bool(v) for k, v in raw.items()}


def _samples_payload(cfg: TxNovaConfig, rows: list[SampleRow], preflight: dict[str, Any]) -> str:
    dups = _dup_map(preflight)
    samples = samples_to_jsonable(rows)
    for s in samples:
        s["dup_flag_seen"] = dups.get(s["sample_id"], False)
    return json.dumps(
        {
            "assembly": cfg.genome.assembly,
            "library_layout": cfg.quantify.library_layout,
            "de_enabled": cfg.de.enabled,
            "treat_min_detected_replicates": cfg.filters.treat_min_detected_replicates,
            "threads": cfg.threads,
            "samples": samples,
        }
    )


def run_pipeline(
    cfg: TxNovaConfig,
    *,
    config_path: Path,
    force: bool = False,
    assembler: Assembler | None = None,
) -> dict[str, Any]:
    rows = load_samples(cfg.samples)
    _apply_contrast_policy(cfg, rows)
    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    remove_orphan_staging(cfg.output_dir)
    if force:
        clear_stamps(cfg.output_dir)
    report = run_preflight(cfg, rows)
    write_preflight_json(cfg, report)
    require_ok(report)

    asm = assembler or Assembler()
    disc = _run_discovery(cfg, rows, report, asm)

    run_json = build_run_json(cfg, rows, report, config_path=config_path, phase="4")
    run_json["assembler"] = {"name": disc["version"], "path": disc["path"]}
    run_json["discovery"] = "phase4"
    warnings = list(run_json.get("warnings") or [])
    warnings.extend(disc.get("warnings") or [])
    run_json["warnings"] = warnings
    run_json["metrics"] = {
        "n_samples": len(rows),
        "n_merged_transcripts": disc["n_merged_transcripts"],
        "n_u": disc["n_u"],
        "n_candidates": disc["n_candidates"],
        "n_structure_error": disc.get("n_structure_error", 0),
        "n_residual_degenerate": disc.get("n_residual_degenerate", 0),
        "dropped_fragments": disc.get("dropped_fragments", 0),
    }
    write_run_json(cfg, run_json)

    n_cand = int(disc["n_candidates"])
    md = render_report(cfg, rows, report, n_candidates=n_cand)
    html = render_html(cfg, rows, report, n_candidates=n_cand)
    write_report(cfg, md, html)
    log.info("run complete: %s candidates in %s", n_cand, cfg.output_dir)
    return report


def _run_discovery(
    cfg: TxNovaConfig,
    rows: list[SampleRow],
    preflight: dict[str, Any],
    asm: Assembler,
) -> dict[str, Any]:
    # Stamp fingerprints may only name files that already exist and that
    # this run will not rewrite before the comparison.
    out = cfg.output_dir
    assembly = out / "assembly"
    classify_dir = out / "classify"
    quant_dir = out / "quantify"
    struct_dir = out / "structure"
    cand_dir = out / "candidates"
    stamps = out / "stamps"
    for d in (assembly, classify_dir, quant_dir, struct_dir, cand_dir, stamps):
        d.mkdir(parents=True, exist_ok=True)

    version = asm.version()
    binary_path = str(asm.resolve_binary())
    core = _core()
    engine = str(core.core_version())
    plan = preflight.get("threads") or {}
    log.info(
        "work plan: BAM ceiling %s (live follows MemAvailable), assembler %s "
        "(cpus=%s usable=%s MemAvailable=%s auto=%s)",
        plan.get("bam_workers"),
        version,
        plan.get("cpus"),
        plan.get("usable"),
        plan.get("available_bytes"),
        plan.get("auto"),
    )

    fp_base = {
        "sheet": sheet_sha(cfg.samples),
        "bams": bam_fingerprint(rows),
        "assembler": version,
        "annotation": str(cfg.genome.annotation),
    }
    assemble_fp = {
        **fp_base,
        "step": "assemble",
        "cfg": {"assembler": ASSEMBLER_VERSION},
    }
    assemble_stamp = stamps / "assemble.json"
    sample_gtfs: list[Path] = []
    if stamp_matches(assemble_stamp, assemble_fp) and all(
        (assembly / f"{r.sample_id}.gtf").is_file() for r in rows
    ):
        log.info("stamp hit assemble")
        sample_gtfs = [assembly / f"{r.sample_id}.gtf" for r in rows]
    else:
        for r in rows:
            dest = assembly / f"{r.sample_id}.gtf"
            asm.assemble(
                r.bam,
                cfg.genome.annotation,
                dest,
                strandedness=r.strandedness,
                threads=1,
                extra_args=[],
            )
            sample_gtfs.append(dest)
        write_stamp(assemble_stamp, assemble_fp)

    merge_fp = {
        **fp_base,
        "step": "merge",
        "annotation": path_fingerprint([cfg.genome.annotation]),
        "cfg": {"assembler": ASSEMBLER_VERSION},
    }
    merge_stamp = stamps / "merge.json"
    merged = assembly / "merged.gtf"
    if not (stamp_matches(merge_stamp, merge_fp) and merged.is_file()):
        asm.merge(sample_gtfs, cfg.genome.annotation, merged, extra_args=[])
        write_stamp(merge_stamp, merge_fp)
    else:
        log.info("stamp hit merge")

    class_tsv = classify_dir / "transcripts.class.tsv"
    reps = classify_dir / "representatives.tsv"
    class_cfg = json.dumps({"strandedness": rows[0].strandedness})

    payload = _samples_payload(cfg, rows, preflight)
    leak_path = cand_dir / LEAK_FILENAME
    leak_scan_fp = {
        **fp_base,
        "step": "leak",
        "engine": engine,
        "inputs": path_fingerprint([merged]),
        "cfg": {
            "min_mapq": cfg.quantify.min_mapq,
            "require_unique_nh": cfg.quantify.require_unique_nh,
            "skip_duplicate": cfg.quantify.skip_duplicate,
            "library_layout": cfg.quantify.library_layout,
            "min_samples": 2,
            "min_support": 2,
            "max_control_support": 0,
            "cohorts": ["silent", "shared", "cohort"],
        },
    }
    if not (stamp_matches(stamps / "leak.json", leak_scan_fp) and leak_path.is_file()):
        leak_cfg = {
            "min_mapq": cfg.quantify.min_mapq,
            "require_unique_nh": cfg.quantify.require_unique_nh,
            "strandedness": rows[0].strandedness,
            "skip_duplicate": cfg.quantify.skip_duplicate,
            "library_layout": cfg.quantify.library_layout,
            "threads": cfg.threads,
            "min_samples": 2,
            "min_support": 2,
            "max_control_support": 0,
        }
        core.leak_scan(
            str(merged),
            payload,
            str(leak_path),
            json.dumps(leak_cfg),
        )
        write_stamp(stamps / "leak.json", leak_scan_fp)
    else:
        log.info("stamp hit leak")

    residual_path = cand_dir / RESIDUAL_FILENAME
    residual_gtf = cand_dir / RESIDUAL_GTF
    residual_gtfs = [cfg.genome.annotation]
    if cfg.genome.naming_annotation is not None:
        residual_gtfs.append(cfg.genome.naming_annotation)
    residual_fp = {
        **fp_base,
        "step": "residual",
        "inputs": path_fingerprint([leak_path, *residual_gtfs]),
        "cfg": {
            "min_exon_nt": MIN_EXON_NT,
            "max_exon_nt": MAX_EXON_NT,
            "min_terminal_nt": MIN_TERMINAL_NT,
            "max_terminal_nt": MAX_TERMINAL_NT,
            "terminal": "cohort_coverage",
            "either_strand_gene_body": True,
            "close_same_strand_nt": CLOSE_SAME_STRAND_NT,
            "min_intron_nt": MIN_INTRON_NT,
            "small_rna_nt": CLOSE_SAME_STRAND_NT,
            "keep_stub_terminal": True,
            "clip_contig": True,
            "no_chain_across_gene": True,
        },
    }
    if not (
        stamp_matches(stamps / "residual.json", residual_fp)
        and residual_path.is_file()
        and residual_gtf.is_file()
        and (cand_dir / RESIDUAL_STATS).is_file()
    ):
        write_residuals(
            leak_path,
            residual_path,
            gene_gtfs=residual_gtfs,
            samples_json=payload,
            min_mapq=cfg.quantify.min_mapq,
            skip_duplicate=cfg.quantify.skip_duplicate,
            fasta=cfg.genome.fasta,
        )
        write_stamp(stamps / "residual.json", residual_fp)
    else:
        log.info("stamp hit residual")

    universe = assembly / UNIVERSE_GTF
    universe_fp = {
        **fp_base,
        "step": "universe",
        "inputs": path_fingerprint([merged, residual_gtf]),
    }
    if not (stamp_matches(stamps / "universe.json", universe_fp) and universe.is_file()):
        write_universe_gtf(merged, residual_gtf, universe)
        write_stamp(stamps / "universe.json", universe_fp)
    else:
        log.info("stamp hit universe")

    classify_fp = {
        **fp_base,
        "step": "classify",
        "engine": engine,
        "merged": path_fingerprint([universe]),
        "cfg": class_cfg,
    }
    if (
        stamp_matches(stamps / "classify.json", classify_fp)
        and class_tsv.is_file()
        and reps.is_file()
    ):
        log.info("stamp hit classify")
    else:
        core.classify_gtfs(
            str(universe),
            str(merged),
            str(class_tsv),
            class_cfg,
        )
        pick_representative_transcript(class_tsv, reps)
        write_stamp(stamps / "classify.json", classify_fp)

    qcfg_obj = {
        "min_mapq": cfg.quantify.min_mapq,
        "require_unique_nh": cfg.quantify.require_unique_nh,
        "strandedness": rows[0].strandedness,
        "skip_duplicate": cfg.quantify.skip_duplicate,
        "library_layout": cfg.quantify.library_layout,
        "threads": cfg.threads,
    }
    qcfg = json.dumps(qcfg_obj)
    quantify_fp = {
        **fp_base,
        "step": "quantify",
        "engine": engine,
        "merged": path_fingerprint([universe]),
        "cfg": qcfg_obj,
    }
    dropped_fragments = 0
    if not (
        stamp_matches(stamps / "quantify.json", quantify_fp)
        and (quant_dir / "locus_tpm.tsv").is_file()
        and (quant_dir / "locus_counts.tsv").is_file()
    ):
        qout = core.quantify_gtf(str(universe), payload, str(quant_dir), qcfg)
        dropped_fragments = int((qout or {}).get("dropped_fragments") or 0)
        write_stamp(
            stamps / "quantify.json",
            quantify_fp,
            outputs={"dropped_fragments": dropped_fragments},
        )
    else:
        log.info("stamp hit quantify")
        dropped_fragments = int(
            stamp_outputs(stamps / "quantify.json").get("dropped_fragments") or 0
        )

    scfg_obj = {
        "strandedness": rows[0].strandedness,
        "min_mapq": cfg.quantify.min_mapq,
        "skip_duplicate": cfg.quantify.skip_duplicate,
        "discontinuity_window_bp": cfg.filters.discontinuity_window_bp,
        "discontinuity_valley_bp": cfg.filters.discontinuity_valley_bp,
        "threads": cfg.threads,
    }
    scfg = json.dumps(scfg_obj)
    struct_tsv = struct_dir / "structure.features.tsv"
    structure_fp = {
        **fp_base,
        "step": "structure",
        "engine": engine,
        "reps": path_fingerprint([reps]),
        "fasta": str(cfg.genome.fasta),
        "cfg": scfg_obj,
    }
    if not (stamp_matches(stamps / "structure.json", structure_fp) and struct_tsv.is_file()):
        core.structure_scan(
            str(universe),
            str(merged),
            str(class_tsv),
            str(reps),
            str(cfg.genome.fasta),
            payload,
            str(struct_tsv),
            scfg,
        )
        write_stamp(stamps / "structure.json", structure_fp)
    else:
        log.info("stamp hit structure")

    gates_path = cand_dir / "candidates.gates.tsv"
    unnamed_path = cand_dir / UNNAMED_FILENAME
    de_view = cand_dir / "candidates.de.tsv"
    cand_path = cand_dir / "candidates.tsv"
    gates_fp = {
        **fp_base,
        "step": "gates",
        "engine": engine,
        "inputs": path_fingerprint(
            [
                class_tsv,
                reps,
                struct_tsv,
                quant_dir / "locus_tpm.tsv",
                quant_dir / "locus_counts.tsv",
            ]
        ),
        "cfg": cfg.filters.model_dump(mode="json", by_alias=True),
    }
    if not (
        stamp_matches(stamps / "gates.json", gates_fp)
        and gates_path.is_file()
        and unnamed_path.is_file()
    ):
        apply_gates(
            cfg,
            rows,
            class_tsv=class_tsv,
            reps_tsv=reps,
            structure_tsv=struct_tsv,
            locus_tpm=quant_dir / "locus_tpm.tsv",
            locus_counts=quant_dir / "locus_counts.tsv",
            out_dir=cand_dir,
            filename="candidates.gates.tsv",
        )
        attach_assembly_evidence(
            gates_path,
            rows,
            assembly,
            unstranded=all(r.strandedness == "unstranded" for r in rows),
        )
        write_stamp(stamps / "gates.json", gates_fp)
    else:
        log.info("stamp hit gates")

    if cfg.de.enabled:
        from txnova.stats import run_deseq

        de_tsv = quant_dir / "de.tsv"
        de_fp = {
            **fp_base,
            "step": "de",
            "inputs": path_fingerprint([quant_dir / "locus_counts.tsv", gates_path]),
            "cfg": cfg.de.model_dump(mode="json"),
            "policy": "wald_or_independent_filter_na",
        }
        if not (
            stamp_matches(stamps / "de.json", de_fp) and de_view.is_file() and de_tsv.is_file()
        ):
            res = run_deseq(cfg, rows, quant_dir / "locus_counts.tsv", de_tsv)
            apply_de_filter(cfg, rows, gates_path, res, de_view)
            write_stamp(stamps / "de.json", de_fp)
        else:
            log.info("stamp hit de")
        coding_src = de_view
    else:
        # Do not rewrite de_view before the coding fingerprint. Funnel still
        # wants the file; copy after coding uses gates_path as its input.
        coding_src = gates_path
        _copy_candidate_view(gates_path, de_view, rows)

    _apply_naming(cfg, rows, [unnamed_path, gates_path, de_view])

    shared_path = cand_dir / SHARED_FILENAME
    write_shared(
        residual=residual_path,
        structure_tables=[unnamed_path, gates_path],
        dest=shared_path,
    )
    _apply_naming(cfg, rows, [shared_path])

    annotate_leak(
        leak_path,
        gates=gates_path,
        unnamed=unnamed_path,
        candidates=de_view,
        shared=shared_path,
    )

    if cfg.coding.enabled:
        hex_path = resolve_hexamer_table(cfg)
        final_reps = cand_dir / "candidates.reps.tsv"
        coding_cfg = {
            "min_orf_aa": cfg.coding.min_orf_aa,
            "hexamer_table": str(hex_path),
            "hexamer_coding_min": cfg.coding.hexamer_coding_min,
            "hexamer_noncoding_max": cfg.coding.hexamer_noncoding_max,
        }
        coding_fp = {
            **fp_base,
            "step": "coding",
            "engine": engine,
            "inputs": path_fingerprint([coding_src, unnamed_path]),
            "cfg": coding_cfg,
        }
        if not (
            stamp_matches(stamps / "coding.json", coding_fp)
            and (cand_dir / "orfs.tsv").is_file()
            and cand_path.is_file()
            and unnamed_path.is_file()
        ):
            _write_union_reps([coding_src, unnamed_path], final_reps)
            core.scan_orfs(
                str(cfg.genome.fasta),
                str(universe),
                str(final_reps),
                str(cand_dir / "orfs.tsv"),
                str(cand_dir / "peptides.fa"),
                json.dumps(coding_cfg),
            )
            _join_coding(cfg, coding_src, cand_dir / "orfs.tsv", cand_path, rows)
            _join_coding(
                cfg,
                unnamed_path,
                cand_dir / "orfs.tsv",
                unnamed_path,
                rows,
                filter_orf=False,
            )
            write_stamp(stamps / "coding.json", coding_fp)
        else:
            log.info("stamp hit coding")
        orfs_tsv = cand_dir / "orfs.tsv"
        if shared_path.is_file() and orfs_tsv.is_file():
            _join_coding(
                cfg,
                shared_path,
                orfs_tsv,
                shared_path,
                rows,
                filter_orf=False,
            )
    else:
        _copy_candidate_view(coding_src, cand_path, rows)

    rank_path = cand_dir / GENE_RANK_FILENAME
    rank_inputs = [p for p in (unnamed_path, cand_path, cand_dir / "peptides.fa") if p.is_file()]
    rank_fp = {
        **fp_base,
        "step": "gene_rank",
        "inputs": path_fingerprint(rank_inputs),
        "cfg": {"report_n": REPORT_N, "assembly": cfg.genome.assembly},
    }
    if not (stamp_matches(stamps / "gene_rank.json", rank_fp) and rank_path.is_file()):
        write_gene_rank(
            [unnamed_path, cand_path],
            rank_path,
            peptides_fa=cand_dir / "peptides.fa" if (cand_dir / "peptides.fa").is_file() else None,
            assembly=cfg.genome.assembly,
        )
        write_stamp(stamps / "gene_rank.json", rank_fp)
    else:
        log.info("stamp hit gene_rank")
    fold_loci = top_rank_ids(rank_path, REPORT_N)

    if cfg.coding.enabled and cfg.coding.fold:
        fold_dir = out / "fold"
        peptides_fa = cand_dir / "peptides.fa"
        fold_fp = {
            **fp_base,
            "step": "fold",
            "inputs": path_fingerprint([peptides_fa, rank_path]),
            "cfg": {"min_orf_aa": cfg.coding.min_orf_aa, "fold": True, "top_n": REPORT_N},
        }
        if not (
            stamp_matches(stamps / "fold.json", fold_fp) and (fold_dir / "index.html").is_file()
        ):
            fold_peptides(
                peptides_fa,
                [unnamed_path, cand_path],
                fold_dir,
                species=cfg.species,
                min_aa=cfg.coding.min_orf_aa,
                only=fold_loci,
            )
            write_stamp(stamps / "fold.json", fold_fp)
        else:
            log.info("stamp hit fold")
        fn_fp = {
            **fp_base,
            "step": "function",
            "inputs": path_fingerprint([fold_dir / "structures.tsv"]),
        }
        if not (
            stamp_matches(stamps / "function.json", fn_fp) and (fold_dir / "function.tsv").is_file()
        ):
            annotate_structures(fold_dir)
            write_stamp(stamps / "function.json", fn_fp)
        else:
            log.info("stamp hit function")
    if cfg.coding.enabled and cfg.coding.orphan:
        orphan_path = cand_dir / "orphan.tsv"
        peptides_fa = cand_dir / "peptides.fa"
        orphan_fp = {
            **fp_base,
            "step": "orphan",
            "inputs": path_fingerprint(
                [p for p in (rank_path, unnamed_path, cand_path, peptides_fa) if p.is_file()]
            ),
            "cfg": {"orphan": True, "assembly": cfg.genome.assembly, "top_n": REPORT_N},
        }
        if not (stamp_matches(stamps / "orphan.json", orphan_fp) and orphan_path.is_file()):
            annotate_orphans(
                [unnamed_path, cand_path],
                peptides_fa if peptides_fa.is_file() else None,
                orphan_path,
                assembly=cfg.genome.assembly,
                min_orf_aa=cfg.coding.min_orf_aa,
                only=fold_loci,
            )
            write_stamp(stamps / "orphan.json", orphan_fp)
        else:
            log.info("stamp hit orphan")

    write_final_transcripts(class_tsv, cand_path, cand_dir / "transcripts.tsv")

    n_cand = 0
    if cand_path.is_file():
        n_cand = max(0, sum(1 for _ in cand_path.open()) - 1)
    n_tx, n_u = _class_metrics(class_tsv)
    n_struct_err = n_structure_error_loci(struct_tsv)
    n_deg = _residual_degenerate(cand_dir / RESIDUAL_STATS)
    warnings: list[str] = []
    if dropped_fragments:
        msg = f"quantify dropped {dropped_fragments} unpaired PE first-mates (not counted)"
        log.warning(msg)
        warnings.append(msg)
    if n_struct_err:
        msg = f"{n_struct_err} representative transcripts skipped (structure_error)"
        log.warning(msg)
        warnings.append(msg)
    if n_deg:
        msg = f"{n_deg} residual loci dropped (n_exons < 2 after terminal clip)"
        log.warning(msg)
        warnings.append(msg)
    return {
        "version": version,
        "path": binary_path,
        "n_merged_transcripts": n_tx,
        "n_u": n_u,
        "n_candidates": n_cand,
        "n_structure_error": n_struct_err,
        "n_residual_degenerate": n_deg,
        "dropped_fragments": dropped_fragments,
        "warnings": warnings,
    }


def _residual_degenerate(path: Path) -> int:
    if not path.is_file():
        return 0
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return 0
    try:
        return int(payload.get("n_degenerate") or 0)
    except (TypeError, ValueError):
        return 0


def _class_metrics(class_tsv: Path) -> tuple[int | None, int | None]:
    if not class_tsv.is_file():
        return None, None
    df = pd.read_csv(class_tsv, sep="\t")
    n_tx = int(len(df))
    n_u = int((df["class"] == "u").sum()) if "class" in df.columns else None
    return n_tx, n_u


def _copy_candidate_view(src: Path, dest: Path, rows: list[SampleRow]) -> None:
    cand = pd.read_csv(src, sep="\t") if src.is_file() else pd.DataFrame()
    write_candidates(cand, dest, rows)


def n_structure_error_loci(struct_tsv: Path) -> int:
    if not struct_tsv.is_file():
        return 0
    df = pd.read_csv(struct_tsv, sep="\t")
    if df.empty or "structure_error" not in df.columns or "locus_id" not in df.columns:
        return 0
    bad = df.loc[df["structure_error"].map(is_structure_error), "locus_id"]
    return int(bad.nunique())


def write_final_transcripts(class_tsv: Path, cand_path: Path, dest: Path) -> None:
    if not class_tsv.is_file():
        dest.write_text("", encoding="utf-8")
        return
    cls = pd.read_csv(class_tsv, sep="\t")
    cand = pd.read_csv(cand_path, sep="\t") if cand_path.is_file() else pd.DataFrame()
    ids = set(cand["locus_id"]) if not cand.empty and "locus_id" in cand.columns else set()
    iso = cls[cls["gene_id"].isin(ids)] if ids and "gene_id" in cls.columns else cls.iloc[0:0]
    dest.parent.mkdir(parents=True, exist_ok=True)
    iso.to_csv(dest, sep="\t", index=False)


def _write_final_reps(cand_path: Path, dest: Path) -> None:
    _write_union_reps([cand_path], dest)


def _write_union_reps(paths: list[Path], dest: Path) -> None:
    frames: list[pd.DataFrame] = []
    for path in paths:
        if not path.is_file():
            continue
        cand = pd.read_csv(path, sep="\t")
        if cand.empty or "representative_transcript_id" not in cand.columns:
            continue
        frames.append(
            cand[["locus_id", "representative_transcript_id"]].rename(
                columns={"representative_transcript_id": "transcript_id"}
            )
        )
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not frames:
        dest.write_text("locus_id\ttranscript_id\n", encoding="utf-8")
        return
    out = pd.concat(frames, ignore_index=True).drop_duplicates("locus_id")
    out.to_csv(dest, sep="\t", index=False)


def _locus_id_list(path: Path) -> list[str]:
    if not path.is_file():
        return []
    df = pd.read_csv(path, sep="\t")
    if df.empty or "locus_id" not in df.columns:
        return []
    return sorted(df["locus_id"].astype(str))


def _apply_naming(cfg: TxNovaConfig, rows: list[SampleRow], paths: list[Path]) -> None:
    gtf = cfg.genome.naming_annotation
    if gtf is None:
        return
    if not gtf.is_file():
        raise TxNovaError(f"genome.naming_annotation not found: {gtf}")
    stamps = cfg.output_dir / "stamps"
    stamps.mkdir(parents=True, exist_ok=True)
    st = gtf.stat()
    fp = {
        "step": "naming",
        "gtf": {"path": str(gtf), "size": st.st_size, "mtime": int(st.st_mtime)},
        "loci": {p.name: _locus_id_list(p) for p in paths},
    }
    have_col = all(
        (not p.is_file()) or ("named_overlap" in pd.read_csv(p, sep="\t", nrows=0).columns)
        for p in paths
    )
    if stamp_matches(stamps / "naming.json", fp) and have_col:
        log.info("stamp hit naming")
        return
    log.info("naming loci against %s", gtf)
    genes = load_gene_bodies(gtf)
    for path in paths:
        if path.is_file():
            label_candidate_file(path, rows, genes)
    write_stamp(stamps / "naming.json", fp)


def apply_de_filter(
    cfg: TxNovaConfig,
    rows: list[SampleRow],
    gates_path: Path,
    res: pd.DataFrame,
    dest: Path,
) -> None:
    from txnova.stats import de_pass_mask, de_status_series

    cand = pd.read_csv(gates_path, sep="\t") if gates_path.is_file() else pd.DataFrame()
    cand = drop_join_columns(cand, DE_COLUMNS)
    if cand.empty:
        write_candidates(cand, dest, rows)
        return
    de_cols = [
        c for c in ("locus_id", "baseMean", "log2FoldChange", "pvalue", "padj") if c in res.columns
    ]
    merged_de = cand.merge(res[de_cols], on="locus_id", how="left")
    merged_de["de_status"] = de_status_series(merged_de, cfg)
    merged_de["de_pass"] = de_pass_mask(merged_de, cfg)
    merged_de = merged_de[merged_de["de_pass"]].copy()
    write_candidates(merged_de, dest, rows)


def _join_coding(
    cfg: TxNovaConfig,
    src_path: Path,
    orfs_path: Path,
    dest_path: Path,
    rows: list[SampleRow],
    *,
    filter_orf: bool = True,
) -> None:
    cand = pd.read_csv(src_path, sep="\t") if src_path.is_file() else pd.DataFrame()
    cand = drop_join_columns(cand, CODING_COLUMNS)
    if cand.empty or not orfs_path.is_file():
        write_candidates(cand, dest_path, rows)
        return
    orfs = pd.read_csv(orfs_path, sep="\t")
    cols = ["locus_id", *CODING_COLUMNS]
    missing = [c for c in cols if c not in orfs.columns]
    if missing:
        raise TxNovaError(f"{orfs_path} missing columns: {missing}")
    out = cand.merge(orfs[cols], on="locus_id", how="left")
    out["coding_label"] = out["coding_label"].fillna("noncoding")
    for col in ("coding_score", "fickett_score"):
        out[col] = out[col].apply(lambda v: "NA" if pd.isna(v) else v)
    if filter_orf and cfg.coding.require_orf:
        aa = pd.to_numeric(out["longest_orf_aa"], errors="coerce").fillna(0)
        complete = out["orf_complete"].astype(str).str.lower().isin(["true", "1"])
        out = out[complete & (aa >= cfg.coding.min_orf_aa)].copy()
    write_candidates(out, dest_path, rows)


def write_report_from_outputs(cfg: TxNovaConfig) -> Path:
    rows = load_samples(cfg.samples)
    pf = cfg.output_dir / "preflight.json"
    if not pf.is_file():
        raise TxNovaError(f"no preflight.json in {cfg.output_dir}; run `txnova run` first")
    report = json.loads(pf.read_text(encoding="utf-8"))
    cand = cfg.output_dir / "candidates" / "candidates.tsv"
    n = 0
    if cand.is_file():
        n = max(0, sum(1 for _ in cand.open()) - 1)
    md = render_report(cfg, rows, report, n_candidates=n)
    html = render_html(cfg, rows, report, n_candidates=n)
    return write_report(cfg, md, html)
