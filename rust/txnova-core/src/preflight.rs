use crate::bam;
use crate::error::CoreError;
use crate::fasta;
use crate::gtf;
use crate::sys_mem;
use crate::types::{PreflightInput, PreflightReport, SampleOut, ThreadPlan};
use serde_json::json;
use std::collections::{BTreeSet, HashMap};
use std::path::Path;

const SORT_SCAN: usize = 10_000;
const LAYOUT_SCAN: usize = 100_000;

pub fn run_preflight(samples_json: &str, fasta: &str, annotation_gtf: &str) -> PreflightReport {
    match run_preflight_inner(samples_json, fasta, annotation_gtf) {
        Ok(r) => r,
        Err((mut r, err)) => {
            r.ok = false;
            if !r.errors.iter().any(|e| e == &err) {
                r.errors.push(err);
            }
            r
        }
    }
}

fn fail_report(
    report: PreflightReport,
    err: String,
) -> std::result::Result<PreflightReport, (PreflightReport, String)> {
    Err((report, err))
}

fn run_preflight_inner(
    samples_json: &str,
    fasta: &str,
    annotation_gtf: &str,
) -> std::result::Result<PreflightReport, (PreflightReport, String)> {
    let mut report = PreflightReport {
        assembly: String::new(),
        ..PreflightReport::default()
    };
    let input_early: std::result::Result<PreflightInput, serde_json::Error> =
        serde_json::from_str(samples_json);
    let requested_threads = input_early.as_ref().map(|i| i.threads).unwrap_or(0);
    let n_samples = input_early.as_ref().map(|i| i.samples.len()).unwrap_or(0);
    let plan = sys_mem::resolve_work_plan(
        requested_threads,
        sys_mem::cpu_count(),
        sys_mem::available_memory_bytes(),
        n_samples,
    );
    report.threads = ThreadPlan {
        requested: plan.requested,
        cpus: plan.cpus,
        usable: plan.usable,
        available_bytes: plan.available_bytes,
        bam_workers: plan.bam_workers,
        auto: plan.auto,
    };

    let input: PreflightInput = match serde_json::from_str(samples_json) {
        Ok(v) => v,
        Err(e) => return fail_report(report, format!("invalid samples_json: {e}")),
    };
    report.assembly = input.assembly.clone();

    if input.samples.is_empty() {
        return fail_report(report, "sample sheet is empty".into());
    }

    let mut groups = HashMap::<String, usize>::new();
    let mut strands = BTreeSet::<String>::new();
    for s in &input.samples {
        if s.sample_id.is_empty() {
            return fail_report(report, "sample_id is empty".into());
        }
        if !s.group.is_empty() && s.group != "control" && s.group != "treat" {
            return fail_report(
                report,
                format!(
                    "sample {} has group {:?}; only control|treat are accepted (no case/ko/wt aliases)",
                    s.sample_id, s.group
                ),
            );
        }
        if s.strandedness != "unstranded" && s.strandedness != "fr" && s.strandedness != "rf" {
            return fail_report(
                report,
                format!(
                    "sample {} has strandedness {:?}; only unstranded|fr|rf",
                    s.sample_id, s.strandedness
                ),
            );
        }
        *groups.entry(s.group.clone()).or_insert(0) += 1;
        strands.insert(s.strandedness.clone());
    }
    let n_control = *groups.get("control").unwrap_or(&0);
    let n_treat = *groups.get("treat").unwrap_or(&0);
    report.n_control = n_control;
    report.n_treat = n_treat;
    if input.samples.len() < 2 {
        return fail_report(
            report,
            format!(
                "need at least 2 samples for residual harvest (got {})",
                input.samples.len()
            ),
        );
    }
    if strands.len() != 1 {
        return fail_report(
            report,
            format!("all samples must share the same strandedness; got {strands:?}"),
        );
    }
    report.strandedness = strands.iter().next().unwrap().clone();
    if input.de_enabled && (n_control < 2 || n_treat < 2) {
        return fail_report(
            report,
            format!(
                "de.enabled requires ≥2 control and ≥2 treat (got control={n_control} treat={n_treat})"
            ),
        );
    }

    let fasta_path = Path::new(fasta);
    let fai = match fasta::load_fai(fasta_path) {
        Ok(v) => v,
        Err(e) => return fail_report(report, e.to_string()),
    };
    let fai_map = fasta::fai_map(&fai);

    let gtf_names = match gtf::gtf_seqnames(Path::new(annotation_gtf)) {
        Ok(v) => v,
        Err(e) => return fail_report(report, e.to_string()),
    };

    let mut common_sq: Option<Vec<(String, u64)>> = None;
    let mut families = BTreeSet::<String>::new();
    let mut layouts = BTreeSet::<String>::new();
    let requested_layout = input.library_layout.as_str();

    for sample in &input.samples {
        let bam_path = Path::new(&sample.bam);
        if !bam_path.exists() {
            return fail_report(
                report,
                format!(
                    "sample {}: BAM not found: {}",
                    sample.sample_id,
                    bam_path.display()
                ),
            );
        }
        if let Err(e) = bam::validate_bam_integrity(bam_path) {
            return fail_report(report, format!("sample {}: {e}", sample.sample_id));
        }
        let ireader = match bam::open_indexed(bam_path) {
            Ok(r) => r,
            Err(e) => return fail_report(report, format!("sample {}: {e}", sample.sample_id)),
        };
        let text = bam::header_text(&ireader);
        drop(ireader);

        let sq = bam::parse_sq(&text);
        if sq.is_empty() {
            return fail_report(
                report,
                format!("sample {}: BAM has no @SQ", sample.sample_id),
            );
        }

        for (sn, ln) in &sq {
            match fai_map.get(sn) {
                None => {
                    let bam_ex: Vec<_> = sq.iter().take(5).map(|(s, _)| s.clone()).collect();
                    let fa_ex: Vec<_> = fai.iter().take(5).map(|(s, _)| s.clone()).collect();
                    return fail_report(
                        report,
                        format!(
                            "sample {}: BAM contig {sn} is not in FASTA .fai \
                             (BAM e.g. {bam_ex:?}, FASTA e.g. {fa_ex:?}). Names must match literally.",
                            sample.sample_id
                        ),
                    );
                }
                Some(fa_ln) if fa_ln != ln => {
                    return fail_report(
                        report,
                        format!(
                            "sample {}: contig {sn} LN={ln} but FASTA .fai has {fa_ln}",
                            sample.sample_id
                        ),
                    );
                }
                _ => {}
            }
        }

        match &common_sq {
            None => common_sq = Some(sq.clone()),
            Some(prev) if !sq_equal(prev, &sq) => {
                return fail_report(
                    report,
                    format!(
                        "sample {}: BAM @SQ (SN+LN) does not match the other samples",
                        sample.sample_id
                    ),
                );
            }
            _ => {}
        }

        let mut seq = match bam::open_sequential(bam_path) {
            Ok(r) => r,
            Err(e) => return fail_report(report, format!("sample {}: {e}", sample.sample_id)),
        };
        if let Err(e) = bam::check_coordinate_order_seq(&mut seq, SORT_SCAN) {
            return fail_report(report, format!("sample {}: {e}", sample.sample_id));
        }

        let family = match bam::infer_aligner_family(&text) {
            Ok(f) => f,
            Err(e) => return fail_report(report, format!("sample {}: {e}", sample.sample_id)),
        };
        families.insert(family.as_str().to_string());

        let mut seq = match bam::open_sequential(bam_path) {
            Ok(r) => r,
            Err(e) => return fail_report(report, format!("sample {}: {e}", sample.sample_id)),
        };
        let layout = match bam::scan_layout_seq(&mut seq, LAYOUT_SCAN) {
            Ok(v) => v,
            Err(e) => return fail_report(report, format!("sample {}: {e}", sample.sample_id)),
        };
        if layout.mixed {
            return fail_report(
                report,
                format!(
                    "sample {}: BAM mixes paired and single-end primary mapped reads (>1% each)",
                    sample.sample_id
                ),
            );
        }
        let observed = if layout.paired { "paired" } else { "single" };
        match requested_layout {
            "auto" => {}
            "paired" | "single" if requested_layout != observed => {
                return fail_report(
                    report,
                    format!(
                        "sample {}: quantify.library_layout={requested_layout} but BAM looks {observed}",
                        sample.sample_id
                    ),
                );
            }
            "paired" | "single" => {}
            other => {
                return fail_report(
                    report,
                    format!("quantify.library_layout must be auto|paired|single, got {other}"),
                );
            }
        }
        layouts.insert(observed.to_string());
        report
            .dup_flag_seen
            .insert(sample.sample_id.clone(), json!(layout.dup_seen));
        if !layout.dup_seen {
            report.warnings.push(format!(
                "sample {}: no 0x400 duplicate flag in the first {LAYOUT_SCAN} primary mapped reads; \
                 quantify.skip_duplicate=auto will count all reads",
                sample.sample_id
            ));
        }
        report.samples.push(SampleOut {
            sample_id: sample.sample_id.clone(),
            n_sq: sq.len(),
        });
    }

    if families.len() != 1 {
        return fail_report(
            report,
            format!(
                "samples mix aligner families {families:?}; all BAMs must be STAR or all HISAT2"
            ),
        );
    }
    if layouts.len() != 1 {
        return fail_report(report, format!("samples mix library layouts {layouts:?}"));
    }
    report.aligner_family = families.iter().next().unwrap().clone();
    report.library_layout = layouts.iter().next().unwrap().clone();

    let sq = common_sq.unwrap();
    report.sq_sha256 = bam::sq_sha256(&sq);
    let bam_names: BTreeSet<String> = sq.iter().map(|(s, _)| s.clone()).collect();
    match gtf::check_gtf_vs_bam(&gtf_names, &bam_names) {
        Ok(n) => report.gtf_seqnames = n,
        Err(e) => return fail_report(report, e.to_string()),
    }

    report.ok = true;
    Ok(report)
}

fn sq_equal(a: &[(String, u64)], b: &[(String, u64)]) -> bool {
    let mut aa = a.to_vec();
    let mut bb = b.to_vec();
    aa.sort();
    bb.sort();
    aa == bb
}

#[allow(dead_code)]
pub fn core_error_from(e: CoreError) -> String {
    e.to_string()
}
