from __future__ import annotations

from pathlib import Path

from txnova.fold import (
    mean_ca_b_factor,
    parse_peptides_fa,
    render_structure_html,
    taxon_for_species,
)


def test_taxon() -> None:
    assert taxon_for_species("mouse") == 10090
    assert taxon_for_species("human") == 9606


def test_parse_peptides(tmp_path: Path) -> None:
    p = tmp_path / "p.fa"
    p.write_text(">MSTRG.1|MSTRG.1.1|3\nAAA\n>MSTRG.2|x|2\nGG\n", encoding="utf-8")
    recs = parse_peptides_fa(p)
    assert recs == {"MSTRG.1": "AAA", "MSTRG.2": "GG"}


def test_mean_ca_b_factor() -> None:
    pdb = (
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 40.00\n"
        "ATOM      2  CA  ALA A   1       1.000   0.000   0.000  1.00 80.00\n"
        "ATOM      3  CA  GLY A   2       2.000   0.000   0.000  1.00 90.00\n"
    )
    assert mean_ca_b_factor(pdb) == 85.0


def test_html_embeds_pdb_and_plddt_palette() -> None:
    html = render_structure_html(
        title="MSTRG.1",
        pdb_text="ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 95.00\n",
        source="esmfold",
        mean_plddt=95.0,
    )
    assert "3Dmol" in html
    assert "#0053D6" in html
    assert "MSTRG.1" in html
    assert "ATOM" in html
