mod bam;
mod classcode;
mod coding;
mod coverage;
mod error;
mod fasta;
mod gtf;
mod interval;
mod leak;
mod orf;
mod preflight;
mod quantify;
mod splice;
mod structure;
mod sys_mem;
mod terminal;
mod types;

pub use preflight::run_preflight;
pub use sys_mem::{resolve_work_plan, WorkPlan};

use pyo3::create_exception;
use pyo3::exceptions::PyException;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyModule};

create_exception!(txnova_core, TxNovaError, PyException);

#[pyfunction]
fn core_version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

fn json_dict(py: Python<'_>, value: serde_json::Value) -> PyResult<Py<PyDict>> {
    let json =
        serde_json::to_string(&value).map_err(|e| PyErr::new::<TxNovaError, _>(e.to_string()))?;
    let json_mod = py.import("json")?;
    json_mod.call_method1("loads", (json,))?.extract()
}

fn map_err(e: error::CoreError) -> PyErr {
    PyErr::new::<TxNovaError, _>(e.to_string())
}

#[pyfunction]
fn preflight_bams(
    py: Python<'_>,
    samples_json: &str,
    fasta: &str,
    annotation_gtf: &str,
) -> PyResult<Py<PyDict>> {
    let report = preflight::run_preflight(samples_json, fasta, annotation_gtf);
    let value =
        serde_json::to_value(&report).map_err(|e| PyErr::new::<TxNovaError, _>(e.to_string()))?;
    json_dict(py, value)
}

#[pyfunction]
fn classify_gtfs(
    py: Python<'_>,
    merged_gtf: &str,
    ref_gtf: &str,
    out_tsv: &str,
    cfg_json: &str,
) -> PyResult<Py<PyDict>> {
    let (n_transcripts, n_u) =
        classcode::classify_gtfs(merged_gtf, ref_gtf, out_tsv, cfg_json).map_err(map_err)?;
    json_dict(
        py,
        serde_json::json!({ "n_transcripts": n_transcripts, "n_u": n_u }),
    )
}

#[pyfunction]
fn structure_scan(
    py: Python<'_>,
    merged_gtf: &str,
    ref_gtf: &str,
    class_tsv: &str,
    representatives_tsv: &str,
    fasta: &str,
    samples_json: &str,
    out_tsv: &str,
    cfg_json: &str,
) -> PyResult<Py<PyDict>> {
    let n = structure::structure_scan(
        merged_gtf,
        ref_gtf,
        class_tsv,
        representatives_tsv,
        fasta,
        samples_json,
        out_tsv,
        cfg_json,
    )
    .map_err(map_err)?;
    json_dict(py, serde_json::json!({ "n_loci_scanned": n }))
}

#[pyfunction]
fn quantify_gtf(
    py: Python<'_>,
    merged_gtf: &str,
    samples_json: &str,
    out_dir: &str,
    cfg_json: &str,
) -> PyResult<Py<PyDict>> {
    let (n_transcripts, n_loci, dropped_fragments) =
        quantify::quantify_gtf(merged_gtf, samples_json, out_dir, cfg_json).map_err(map_err)?;
    json_dict(
        py,
        serde_json::json!({
            "n_transcripts": n_transcripts,
            "n_loci": n_loci,
            "dropped_fragments": dropped_fragments,
        }),
    )
}

#[pyfunction]
fn leak_scan(
    py: Python<'_>,
    merged_gtf: &str,
    samples_json: &str,
    out_tsv: &str,
    cfg_json: &str,
) -> PyResult<Py<PyDict>> {
    let (n, dropped_fragments) =
        leak::leak_scan(merged_gtf, samples_json, out_tsv, cfg_json).map_err(map_err)?;
    json_dict(
        py,
        serde_json::json!({ "n_leak": n, "dropped_fragments": dropped_fragments }),
    )
}

#[pyfunction]
fn scan_orfs(
    py: Python<'_>,
    fasta: &str,
    merged_gtf: &str,
    representatives_tsv: &str,
    out_tsv: &str,
    peptides_fa: &str,
    cfg_json: &str,
) -> PyResult<Py<PyDict>> {
    let n = orf::scan_orfs(
        fasta,
        merged_gtf,
        representatives_tsv,
        out_tsv,
        peptides_fa,
        cfg_json,
    )
    .map_err(map_err)?;
    json_dict(py, serde_json::json!({ "n_orfs": n }))
}

#[pyfunction]
fn extend_terminals(
    py: Python<'_>,
    loci_tsv: &str,
    out_tsv: &str,
    samples_json: &str,
    cfg_json: &str,
) -> PyResult<Py<PyDict>> {
    let n =
        terminal::extend_terminals(loci_tsv, out_tsv, samples_json, cfg_json).map_err(map_err)?;
    json_dict(py, serde_json::json!({ "n_loci": n }))
}

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("TxNovaError", m.py().get_type::<TxNovaError>())?;
    m.add_function(wrap_pyfunction!(core_version, m)?)?;
    m.add_function(wrap_pyfunction!(preflight_bams, m)?)?;
    m.add_function(wrap_pyfunction!(classify_gtfs, m)?)?;
    m.add_function(wrap_pyfunction!(structure_scan, m)?)?;
    m.add_function(wrap_pyfunction!(quantify_gtf, m)?)?;
    m.add_function(wrap_pyfunction!(scan_orfs, m)?)?;
    m.add_function(wrap_pyfunction!(leak_scan, m)?)?;
    m.add_function(wrap_pyfunction!(extend_terminals, m)?)?;
    Ok(())
}
