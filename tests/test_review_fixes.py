from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest
from txnova.config import TxNovaConfig, load_config, peek_gtf_species
from txnova.errors import TxNovaError
from txnova.fold import render_structure_embed, render_structure_html
from txnova.gates import apply_gates
from txnova.naming import chrom_key
from txnova.residual import _cluster_record, write_residuals
from txnova.samples import SampleRow


def _member(
    *,
    start: int,
    end: int,
    n_detected: int = 0,
    support_sum: int = 0,
    treat_n_detected: int = 0,
    treat_sum: int = 0,
    control_max: int = 0,
    control_n_detected: int = 0,
    cohort: str = "",
) -> dict:
    return {
        "chrom": "chr8",
        "start": start,
        "end": end,
        "strand": "+",
        "n_detected": n_detected,
        "support_sum": support_sum,
        "treat_n_detected": treat_n_detected,
        "treat_sum": treat_sum,
        "control_max": control_max,
        "control_n_detected": control_n_detected,
        "cohort": cohort,
        "nearest_distance_bp": 20000,
        "nearest_gene_name": "Far",
    }


def test_cluster_record_catalog_only_uses_n_detected_and_keeps_cohort() -> None:
    members = [
        _member(start=1000, end=1100, n_detected=2, support_sum=10, cohort="cohort"),
        _member(start=1200, end=1300, n_detected=3, support_sum=20, cohort="cohort"),
    ]
    rec = _cluster_record(members, "RSDL.1")
    assert rec["treat_n_detected"] == 3
    assert rec["treat_sum"] == 30
    assert rec["cohort"] == "cohort"


def test_cluster_record_contrast_uses_treat_fields_not_union() -> None:
    members = [
        _member(
            start=1000,
            end=1100,
            n_detected=8,
            support_sum=50,
            treat_n_detected=4,
            treat_sum=40,
            control_max=5,
            control_n_detected=4,
            cohort="shared",
        ),
        _member(
            start=1200,
            end=1300,
            n_detected=8,
            support_sum=50,
            treat_n_detected=3,
            treat_sum=20,
            control_max=2,
            control_n_detected=3,
            cohort="shared",
        ),
    ]
    rec = _cluster_record(members, "RSDL.1")
    assert rec["treat_n_detected"] == 4
    assert rec["treat_sum"] == 60
    assert rec["cohort"] == "shared"
    assert rec["control_n_detected"] == 4


def test_cluster_record_silent_keeps_cohort() -> None:
    members = [
        _member(
            start=1000,
            end=1100,
            n_detected=4,
            support_sum=12,
            treat_n_detected=4,
            treat_sum=12,
            cohort="silent",
        )
    ]
    rec = _cluster_record(members, "RSDL.1")
    assert rec["treat_n_detected"] == 4
    assert rec["treat_sum"] == 12
    assert rec["cohort"] == "silent"


def test_chrom_key_case_insensitive_chr_prefix() -> None:
    assert chrom_key("CHR8") == "8"
    assert chrom_key("Chr8") == "8"
    assert chrom_key("chr8") == "8"
    assert chrom_key("CHRM") == "MT"
    assert chrom_key("ChrM") == "MT"
    assert chrom_key("chrMT") == "MT"


def test_peek_gtf_species_from_transcript_lines(tmp_path: Path) -> None:
    human = tmp_path / "human_tx.gtf"
    human.write_text(
        'chr1\tHAVANA\ttranscript\t1\t10\t.\t+\t.\tgene_id "ENSG00000139618"; '
        'transcript_id "ENST00000356607";\n',
        encoding="utf-8",
    )
    mouse = tmp_path / "mouse_tx.gtf"
    mouse.write_text(
        'chr1\tHAVANA\ttranscript\t1\t10\t.\t+\t.\tgene_id "ENSMUSG00000000001"; '
        'transcript_id "ENSMUST00000000001";\n',
        encoding="utf-8",
    )
    assert peek_gtf_species(human) == "human"
    assert peek_gtf_species(mouse) == "mouse"


def test_auto_species_raises_when_uninferable(tmp_path: Path) -> None:
    gtf = tmp_path / "unknown.gtf"
    gtf.write_text(
        'chr1\tX\texon\t1\t10\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1";\n',
        encoding="utf-8",
    )
    sheet = tmp_path / "s.tsv"
    sheet.write_text("sample_id\tbam\tgroup\tstrandedness\n", encoding="utf-8")
    cfg_path = tmp_path / "c.yaml"
    cfg_path.write_text(
        f"species: auto\ngenome:\n  fasta: /a.fa\n  annotation: {gtf}\n"
        f"  assembly: UNKNOWN\nsamples: s.tsv\n",
        encoding="utf-8",
    )
    with pytest.raises(TxNovaError, match="could not be inferred"):
        load_config(cfg_path)


def test_explicit_species_still_works_when_gtf_is_uninformative(tmp_path: Path) -> None:
    gtf = tmp_path / "unknown.gtf"
    gtf.write_text(
        'chr1\tX\texon\t1\t10\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1";\n',
        encoding="utf-8",
    )
    sheet = tmp_path / "s.tsv"
    sheet.write_text("sample_id\tbam\tgroup\tstrandedness\n", encoding="utf-8")
    cfg_path = tmp_path / "c.yaml"
    cfg_path.write_text(
        f"species: human\ngenome:\n  fasta: /a.fa\n  annotation: {gtf}\n"
        f"  assembly: UNKNOWN\nsamples: s.tsv\n",
        encoding="utf-8",
    )
    cfg = load_config(cfg_path)
    assert cfg.species == "human"
    assert cfg.genome.assembly == "GRCh38"


def test_discovery_only_gates_passed_omits_detect(tmp_path: Path) -> None:
    class_tsv = tmp_path / "c.tsv"
    class_tsv.write_text(
        "transcript_id\tgene_id\tchrom\tstart\tend\tstrand\tn_exons\tlength_nt\t"
        "exon_structure\tclass\tgene_name\n"
        "T1\tL1\tchr1\t1000\t2000\t+\t2\t400\t1000-1200,1800-2000\tu\t\n",
        encoding="utf-8",
    )
    reps = tmp_path / "reps.tsv"
    from txnova.gates import pick_representative_transcript

    pick_representative_transcript(class_tsv, reps)
    struct = tmp_path / "s.tsv"
    struct.write_text(
        "locus_id\tsample_id\ttranscript_id\tchrom\tstart\tend\tstrand\tn_exons\t"
        "length_nt\tdonors\tacceptors\tcanonical_splice_fraction\tn_introns\t"
        "has_nearest\tnearest_gene_id\tnearest_gene_name\tnearest_distance_bp\t"
        "nearest_strand\tnearest_any_gene_id\tnearest_any_gene_name\t"
        "nearest_any_distance_bp\tnearest_any_strand\tlocus_mean_depth\t"
        "valley_mean\tvalley_possible\tgap_mean_depth\tn_dup_flag_seen\t"
        "junction_support\tbridge_read_count\tstructure_error\n"
        "L1\ta\tT1\tchr1\t1000\t2000\t+\t2\t400\tGT\tAG\t1.0\t1\tfalse\t"
        "\t\t\t\t\t\t\t\t10\t\t\t\t0\t12\t0\n"
        "L1\tb\tT1\tchr1\t1000\t2000\t+\t2\t400\tGT\tAG\t1.0\t1\tfalse\t"
        "\t\t\t\t\t\t\t\t10\t\t\t\t0\t8\t0\n",
        encoding="utf-8",
    )
    tpm = tmp_path / "tpm.tsv"
    tpm.write_text("locus_id\ta\tb\nL1\t5.0\t6.0\n")
    counts = tmp_path / "cnt.tsv"
    counts.write_text("locus_id\ta\tb\nL1\t10\t12\n")
    cfg = TxNovaConfig.model_validate(
        {
            "genome": {"fasta": "/a.fa", "annotation": "/a.gtf"},
            "samples": "/s.tsv",
        }
    )
    rows = [
        SampleRow(sample_id="a", bam=Path("/a.bam"), group="", strandedness="rf", replicate=1),
        SampleRow(sample_id="b", bam=Path("/b.bam"), group="", strandedness="rf", replicate=2),
    ]
    cand = apply_gates(
        cfg,
        rows,
        class_tsv=class_tsv,
        reps_tsv=reps,
        structure_tsv=struct,
        locus_tpm=tpm,
        locus_counts=counts,
        out_dir=tmp_path / "cand",
    )
    assert len(cand) == 1
    passed = str(cand.iloc[0]["gates_passed"])
    assert "detect" not in passed.split(",")
    assert "class" in passed


def test_require_unique_nh_passed_into_extend_cfg(tmp_path: Path, monkeypatch) -> None:
    captured: dict[str, object] = {}

    def fake_extend(loci: str, dest: str, samples_json: str, cfg_json: str) -> None:
        captured["cfg"] = json.loads(cfg_json)
        Path(dest).write_text("", encoding="utf-8")

    import txnova._core as core

    monkeypatch.setattr(core, "extend_terminals", fake_extend)
    src = tmp_path / "leak.tsv"
    pd.DataFrame(
        [
            {
                "chrom": "chr4",
                "start": 1000,
                "end": 1100,
                "strand": "+",
                "status": "unassembled",
                "overlaps_gene": False,
                "treat_sum": 8,
                "treat_n_detected": 3,
                "control_max": 0,
                "nearest_distance_bp": 12000,
                "nearest_gene_name": "Far",
            },
            {
                "chrom": "chr4",
                "start": 1200,
                "end": 1300,
                "strand": "+",
                "status": "unassembled",
                "overlaps_gene": False,
                "treat_sum": 6,
                "treat_n_detected": 3,
                "control_max": 0,
                "nearest_distance_bp": 11800,
                "nearest_gene_name": "Far",
            },
        ]
    ).to_csv(src, sep="\t", index=False)
    dest = tmp_path / "residual.tsv"
    write_residuals(src, dest, samples_json="{}", require_unique_nh=True)
    assert captured["cfg"]["require_unique_nh"] is True
    write_residuals(src, dest, samples_json="{}")
    assert captured["cfg"]["require_unique_nh"] is False


def test_fold_captions_html_escape_title_and_note() -> None:
    title = 'RSDL.1 <script>alert("x")</script>'
    note = "a <b> note & more"
    html = render_structure_embed(
        title=title,
        pdb_text="ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 95.00\n",
        source="esmfold",
        mean_plddt=95.0,
        note=note,
    )
    assert "<script>alert" not in html
    assert "RSDL.1 &lt;script&gt;" in html
    assert "&lt;b&gt;" in html
    assert "&amp;" in html
    page = render_structure_html(
        title=title,
        pdb_text="ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 95.00\n",
        source="esmfold",
        mean_plddt=95.0,
        note=note,
    )
    assert "<script>alert" not in page
    assert "RSDL.1 &lt;script&gt;" in page
