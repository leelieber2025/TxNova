"""Function hints from structure: Foldseek against PDB / AlphaFold DB.

This is fold similarity, not docking. A hit is a homolog with a similar
3D fold; the target description is the function guess. Named proteins
also get the Swiss-Prot FUNCTION comment when a UniProt accession exists.
"""

from __future__ import annotations

import json
import time
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd

from txnova import USER_AGENT
from txnova.logging import get_logger

log = get_logger("txnova.function")

FS_TICKET = "https://search.foldseek.com/api/ticket"
FS_STATUS = "https://search.foldseek.com/api/ticket/{tid}"
FS_RESULT = "https://search.foldseek.com/api/result/{tid}/0"
UNIPROT = "https://rest.uniprot.org/uniprotkb/{acc}.json"

# Keep only confident structural matches
MIN_PROB = 0.5
MAX_EVALUE = 1e-3
TOP_N = 5


def _http_json(url: str, timeout: int = 60) -> dict:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode())


def uniprot_function(acc: str) -> str:
    if not acc:
        return ""
    try:
        data = _http_json(UNIPROT.format(acc=acc), timeout=45)
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError, OSError) as e:
        log.warning("UniProt function failed for %s: %s", acc, e)
        return ""
    bits = []
    for c in data.get("comments") or []:
        if c.get("commentType") != "FUNCTION":
            continue
        for t in c.get("texts") or []:
            v = (t.get("value") or "").strip()
            if v:
                bits.append(v)
    return " ".join(bits)


def _parse_hits(payload: dict) -> list[dict]:
    out: list[dict] = []
    for block in payload.get("results") or []:
        db = str(block.get("db") or "")
        aln = block.get("alignments") or []
        if aln and isinstance(aln[0], list):
            aln = aln[0]
        for h in aln:
            try:
                ev = float(h.get("eval"))
                prob = float(h.get("prob"))
            except (TypeError, ValueError):
                continue
            if prob < MIN_PROB or ev > MAX_EVALUE:
                continue
            target = str(h.get("target") or "")
            desc = target
            # "AF-... ModelName" or "pdbid ... description"
            parts = target.split(None, 1)
            if len(parts) == 2:
                desc = parts[1]
            out.append(
                {
                    "db": db,
                    "target": target.split()[0] if target else "",
                    "description": desc,
                    "prob": prob,
                    "evalue": ev,
                    "seq_id": h.get("seqId"),
                    "taxon": h.get("taxName") or "",
                }
            )
    out.sort(key=lambda r: (r["evalue"], -r["prob"]))
    return out


def _multipart(fields: list[tuple[str, str]], files: list[tuple[str, Path]]) -> tuple[bytes, str]:
    import uuid

    boundary = "----TxNova" + uuid.uuid4().hex
    chunks: list[bytes] = []
    for name, value in fields:
        chunks.append(
            f'--{boundary}\r\nContent-Disposition: form-data; name="{name}"\r\n\r\n{value}\r\n'.encode()
        )
    for name, path in files:
        header = (
            f'--{boundary}\r\nContent-Disposition: form-data; name="{name}"; '
            f'filename="{path.name}"\r\nContent-Type: application/octet-stream\r\n\r\n'
        )
        chunks.append(header.encode())
        chunks.append(path.read_bytes())
        chunks.append(b"\r\n")
    chunks.append(f"--{boundary}--\r\n".encode())
    return b"".join(chunks), boundary


def foldseek_search(pdb: Path, *, timeout_s: int = 180) -> list[dict]:
    body, boundary = _multipart(
        [
            ("mode", "3diaa"),
            ("database[]", "pdb100"),
            ("database[]", "afdb50"),
        ],
        [("q", pdb)],
    )
    req = urllib.request.Request(
        FS_TICKET,
        data=body,
        method="POST",
        headers={
            "User-Agent": USER_AGENT,
            "Content-Type": f"multipart/form-data; boundary={boundary}",
        },
    )
    with urllib.request.urlopen(req, timeout=60) as resp:
        ticket = json.loads(resp.read().decode())
    tid = ticket.get("id")
    if not tid:
        raise RuntimeError(f"no Foldseek ticket: {ticket}")
    deadline = time.time() + timeout_s
    while time.time() < deadline:
        st = _http_json(FS_STATUS.format(tid=tid), timeout=30)
        status = str(st.get("status") or "").upper()
        if status in {"COMPLETE", "COMPLETED"}:
            break
        if status in {"ERROR", "FAILED"}:
            raise RuntimeError(f"Foldseek failed: {st}")
        time.sleep(3)
    else:
        raise TimeoutError(f"Foldseek timed out for {pdb.name}")
    payload = _http_json(FS_RESULT.format(tid=tid), timeout=90)
    return _parse_hits(payload)


def annotate_structures(fold_dir: Path) -> pd.DataFrame:
    table = fold_dir / "structures.tsv"
    if not table.is_file() or table.stat().st_size == 0:
        return pd.DataFrame()
    src = pd.read_csv(table, sep="\t")
    if src.empty:
        return src
    rows: list[dict] = []
    for rec in src.to_dict(orient="records"):
        pdb = Path(str(rec.get("pdb") or ""))
        locus = str(rec.get("locus_id") or "")
        acc = str(rec.get("uniprot") or "")
        if acc in {"", "nan", "NA"}:
            acc = ""
        curated = uniprot_function(acc) if acc else ""
        hits: list[dict] = []
        err = ""
        if pdb.is_file():
            try:
                hits = foldseek_search(pdb)[:TOP_N]
                log.info("%s Foldseek %s hits", locus, len(hits))
            except Exception as e:
                err = str(e)
                log.warning("Foldseek failed for %s: %s", locus, e)
        top = hits[0] if hits else {}
        rows.append(
            {
                "locus_id": locus,
                "gene": rec.get("gene") or "",
                "uniprot": acc,
                "source": rec.get("source") or "",
                "mean_plddt": rec.get("mean_plddt") or "",
                "curated_function": curated,
                "fold_description": top.get("description") or "",
                "fold_target": top.get("target") or "",
                "fold_db": top.get("db") or "",
                "fold_evalue": top.get("evalue") if top else "",
                "fold_prob": top.get("prob") if top else "",
                "fold_seq_id": top.get("seq_id") if top else "",
                "n_hits": len(hits),
                "error": err,
            }
        )
        if hits:
            hit_path = fold_dir / f"{locus}.foldseek.tsv"
            pd.DataFrame(hits).to_csv(hit_path, sep="\t", index=False)
    out = pd.DataFrame(rows)
    out.to_csv(fold_dir / "function.tsv", sep="\t", index=False)
    return out
