"""Protein 3D models for ORFs — AlphaFold DB if named, ESMFold if not.

Visualization is a standalone HTML page with 3Dmol.js, coloured by the
AlphaFold pLDDT bands (same palette as omicverse.mol.view). No docking.
"""

from __future__ import annotations

import json
import re
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Iterable

from txnova.errors import TxNovaError
from txnova.logging import get_logger

log = get_logger("txnova.fold")

MOUSE_TAXON = 10090
HUMAN_TAXON = 9606
ESMFOLD_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"
AF_API = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"
UNIPROT_SEARCH = (
    "https://rest.uniprot.org/uniprotkb/search"
    "?query={query}&fields=accession,gene_primary&format=json&size=1"
)
ESMFOLD_MAX_AA = 400
_UNIPROT_RE = re.compile(
    r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$"
)

# Official AlphaFold pLDDT bands
PLDDT_BANDS = (
    (90.0, 1e9, "#0053D6", "very high (>90)"),
    (70.0, 90.0, "#65CBF3", "high (70–90)"),
    (50.0, 70.0, "#FFDB13", "low (50–70)"),
    (-1e9, 50.0, "#FF7D45", "very low (<50)"),
)


def taxon_for_species(species: str) -> int:
    s = (species or "").strip().lower()
    if s in {"human", "homo sapiens", "hsapiens"}:
        return HUMAN_TAXON
    return MOUSE_TAXON


def parse_peptides_fa(path: Path) -> dict[str, str]:
    """`>locus|transcript|aa` headers from scan_orfs."""
    if not path.is_file():
        return {}
    out: dict[str, str] = {}
    locus = ""
    seq: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if locus and seq:
                out[locus] = "".join(seq)
            locus = line[1:].split("|", 1)[0].strip()
            seq = []
        else:
            seq.append(line.strip())
    if locus and seq:
        out[locus] = "".join(seq)
    return out


def mean_ca_b_factor(pdb_text: str) -> float | None:
    vals: list[float] = []
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM"):
            continue
        if line[12:16].strip() != "CA":
            continue
        try:
            vals.append(float(line[60:66]))
        except ValueError:
            continue
    if not vals:
        return None
    mean = sum(vals) / len(vals)
    # ESMFold sometimes writes 0–1
    if max(vals) <= 1.0:
        mean *= 100.0
    return mean


def _http_get(url: str, timeout: int = 60) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "txnova/0.4"})
    last: Exception | None = None
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, OSError) as e:
            last = e
            if attempt == 2:
                raise
    raise last  # pragma: no cover


def _http_post(url: str, data: bytes, timeout: int = 300) -> bytes:
    last: Exception | None = None
    for attempt in range(3):
        req = urllib.request.Request(
            url,
            data=data,
            method="POST",
            headers={"User-Agent": "txnova/0.4", "Content-Type": "text/plain"},
        )
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, OSError) as e:
            last = e
            if attempt == 2:
                raise
    raise last  # pragma: no cover


def resolve_uniprot(gene: str, taxon: int) -> str | None:
    gene = (gene or "").strip()
    if not gene or gene in {"NA", "nan"}:
        return None
    if _UNIPROT_RE.match(gene):
        return gene
    q = urllib.parse.quote(f"gene_exact:{gene} AND organism_id:{taxon} AND reviewed:true")
    try:
        raw = _http_get(UNIPROT_SEARCH.format(query=q))
        hits = json.loads(raw.decode()).get("results") or []
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError, OSError) as e:
        log.warning("UniProt lookup failed for %s: %s", gene, e)
        return None
    if not hits:
        return None
    return hits[0].get("primaryAccession")


def fetch_alphafold(uniprot: str, dest: Path) -> str:
    recs = json.loads(_http_get(AF_API.format(acc=uniprot)).decode())
    if not recs:
        raise TxNovaError(f"no AlphaFold DB model for {uniprot}")
    url = recs[0].get("pdbUrl")
    if not url:
        raise TxNovaError(f"AlphaFold record for {uniprot} has no pdbUrl")
    dest.write_bytes(_http_get(url))
    return dest.read_text(encoding="utf-8", errors="replace")


def predict_esmfold(sequence: str, dest: Path) -> str:
    seq = "".join(sequence.split()).upper()
    dest.write_bytes(_http_post(ESMFOLD_URL, seq.encode("ascii"), timeout=300))
    return dest.read_text(encoding="utf-8", errors="replace")


def render_structure_html(
    *,
    title: str,
    pdb_text: str,
    source: str,
    mean_plddt: float | None,
    note: str = "",
) -> str:
    pdb_js = json.dumps(pdb_text)
    plddt = "NA" if mean_plddt is None else f"{mean_plddt:.1f}"
    legend = " · ".join(f'<span style="color:{c}">{lab}</span>' for _, _, c, lab in PLDDT_BANDS)
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{title}</title>
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<style>
body {{ margin: 0; font: 15px/1.45 system-ui, sans-serif; color: #1a1a1a; }}
header {{ padding: 0.8rem 1.2rem; border-bottom: 1px solid #ddd; }}
#view {{ width: 100%; height: 640px; }}
.note {{ color: #555; font-size: 0.92rem; }}
</style>
</head>
<body>
<header>
  <strong>{title}</strong>
  <div class="note">source {source} · mean pLDDT {plddt}. {note}</div>
  <div class="note">pLDDT: {legend}</div>
</header>
<div id="view"></div>
<script>
const pdb = {pdb_js};
const el = document.getElementById("view");
const viewer = $3Dmol.createViewer(el, {{backgroundColor: "white"}});
viewer.addModel(pdb, "pdb");
viewer.setStyle({{}}, {{cartoon: {{colorfunc: function(atom) {{
  let b = atom.b;
  if (b <= 1.0) b = b * 100.0;
  if (b >= 90) return "#0053D6";
  if (b >= 70) return "#65CBF3";
  if (b >= 50) return "#FFDB13";
  return "#FF7D45";
}}}}}});
viewer.zoomTo();
viewer.render();
</script>
</body>
</html>
"""


def _write_index(out_dir: Path, rows: list[dict]) -> None:
    lines = [
        "<!DOCTYPE html><html><head><meta charset='utf-8'><title>TxNova structures</title>",
        "<style>body{font:15px/1.45 system-ui,sans-serif;margin:1.5rem} "
        "table{border-collapse:collapse} td,th{border-bottom:1px solid #ddd;"
        "padding:0.3rem 0.6rem;text-align:left}</style></head><body>",
        "<h1>Protein structures</h1>",
        "<p>Named loci use AlphaFold DB. Unnamed ORFs use the ESMFold API. "
        "Colour is pLDDT (blue = confident). Not a functional claim.</p>",
        "<table><tr><th>locus</th><th>source</th><th>uniprot / gene</th>"
        "<th>aa</th><th>mean pLDDT</th><th>view</th></tr>",
    ]
    for r in rows:
        view = r.get("html") or ""
        if not isinstance(view, str) or not view or view == "nan":
            link = ""
        else:
            link = f"<a href='{Path(view).name}'>open</a>"
        lines.append(
            "<tr>"
            f"<td>{r.get('locus_id', '')}</td>"
            f"<td>{r.get('source', '')}</td>"
            f"<td>{r.get('uniprot') or r.get('gene') or ''}</td>"
            f"<td>{r.get('n_aa', '')}</td>"
            f"<td>{r.get('mean_plddt') if r.get('mean_plddt') is not None else ''}</td>"
            f"<td>{link}</td></tr>"
        )
    lines.append("</table></body></html>")
    (out_dir / "index.html").write_text("\n".join(lines), encoding="utf-8")


def fold_peptides(
    peptides_fa: Path,
    tables: Iterable[Path],
    out_dir: Path,
    *,
    species: str,
    min_aa: int,
) -> list[dict]:
    """Write PDB + HTML for each eligible ORF. Network failures are recorded, not fatal."""
    peptides = parse_peptides_fa(peptides_fa)
    meta: dict[str, dict] = {}
    for table in tables:
        if not table.is_file():
            continue
        import pandas as pd

        df = pd.read_csv(table, sep="\t")
        if df.empty or "locus_id" not in df.columns:
            continue
        for rec in df.to_dict(orient="records"):
            meta[str(rec["locus_id"])] = rec

    taxon = taxon_for_species(species)
    out_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []

    for locus, seq in peptides.items():
        info = meta.get(locus, {})
        n_aa = len(seq)
        if n_aa < min_aa:
            continue
        named = str(info.get("named_overlap") or "")
        gene = str(info.get("named_gene_name") or "")
        if gene in {"", "nan", "NA"}:
            gene = ""
        want_af = named == "protein_coding" and bool(gene)
        # Skip loci already named as lncRNA/pseudogene — fold the holes
        # (AlphaFold) and the still-unnamed ORFs (ESMFold).
        if not want_af and named not in {"", "none", "nan"}:
            continue
        source = "alphafold" if want_af else "esmfold"
        if source == "esmfold" and n_aa > ESMFOLD_MAX_AA:
            rows.append(
                {
                    "locus_id": locus,
                    "source": "skipped",
                    "gene": gene,
                    "uniprot": "",
                    "n_aa": n_aa,
                    "mean_plddt": "",
                    "pdb": "",
                    "html": "",
                    "error": f"peptide {n_aa} aa exceeds ESMFold public limit {ESMFOLD_MAX_AA}",
                }
            )
            continue

        pdb_path = out_dir / f"{locus}.pdb"
        html_path = out_dir / f"{locus}.html"
        uniprot = ""
        err = ""
        pdb_text = ""
        try:
            if source == "alphafold":
                uniprot = resolve_uniprot(gene, taxon) or ""
                if not uniprot:
                    source = "esmfold"
                else:
                    pdb_text = fetch_alphafold(uniprot, pdb_path)
            if source == "esmfold":
                pdb_text = predict_esmfold(seq, pdb_path)
        except Exception as e:
            err = str(e)
            log.warning("structure failed for %s: %s", locus, e)
            rows.append(
                {
                    "locus_id": locus,
                    "source": source,
                    "gene": gene,
                    "uniprot": uniprot,
                    "n_aa": n_aa,
                    "mean_plddt": "",
                    "pdb": "",
                    "html": "",
                    "error": err,
                }
            )
            continue

        mean_p = mean_ca_b_factor(pdb_text)
        note = (
            f"AlphaFold DB model for {gene} ({uniprot})."
            if source == "alphafold"
            else "ESMFold prediction of the TxNova ORF. Not a crystal structure."
        )
        html_path.write_text(
            render_structure_html(
                title=locus,
                pdb_text=pdb_text,
                source=source,
                mean_plddt=mean_p,
                note=note,
            ),
            encoding="utf-8",
        )
        rows.append(
            {
                "locus_id": locus,
                "source": source,
                "gene": gene,
                "uniprot": uniprot,
                "n_aa": n_aa,
                "mean_plddt": "" if mean_p is None else round(mean_p, 2),
                "pdb": str(pdb_path),
                "html": str(html_path),
                "error": "",
            }
        )
        log.info("structure %s %s pLDDT=%s", locus, source, rows[-1]["mean_plddt"])

    import pandas as pd

    pd.DataFrame(rows).to_csv(out_dir / "structures.tsv", sep="\t", index=False)
    _write_index(out_dir, rows)
    return rows
