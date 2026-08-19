from __future__ import annotations

from pathlib import Path

import pandas as pd

from txnova.config import TxNovaConfig
from txnova.errors import TxNovaError
from txnova.logging import get_logger
from txnova.samples import SampleRow

log = get_logger("txnova.stats")


def run_deseq(
    cfg: TxNovaConfig,
    rows: list[SampleRow],
    locus_counts: Path,
    out_tsv: Path,
) -> pd.DataFrame:
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except ImportError as e:
        raise TxNovaError(
            "pydeseq2 is a required dependency. Reinstall with `pip install 'txnova'`."
        ) from e

    counts = pd.read_csv(locus_counts, sep="\t").set_index("locus_id")
    sample_ids = [r.sample_id for r in rows]
    missing = [s for s in sample_ids if s not in counts.columns]
    if missing:
        raise TxNovaError(f"locus_counts missing samples {missing}")
    counts = counts[sample_ids].apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)
    keep = counts.sum(axis=1) > 0
    fit = counts.loc[keep]
    empty = counts.loc[~keep]
    if fit.empty:
        res = pd.DataFrame(
            {
                "locus_id": counts.index,
                "baseMean": pd.NA,
                "log2FoldChange": pd.NA,
                "lfcSE": pd.NA,
                "stat": pd.NA,
                "pvalue": pd.NA,
                "padj": pd.NA,
            }
        )
        out_tsv.parent.mkdir(parents=True, exist_ok=True)
        res.to_csv(out_tsv, sep="\t", index=False)
        return res
    # pydeseq2 wants samples x genes
    counts_t = fit.T
    meta = pd.DataFrame(
        {"condition": [r.group for r in rows]},
        index=sample_ids,
    )
    dds = DeseqDataSet(
        counts=counts_t,
        metadata=meta,
        design="~condition",
        quiet=True,
    )
    dds.deseq2()
    stat = DeseqStats(dds, contrast=["condition", "treat", "control"], quiet=True)
    stat.summary()
    res = stat.results_df.copy()
    res.index.name = "locus_id"
    res = res.reset_index()
    for col in ("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"):
        if col not in res.columns:
            res[col] = pd.NA
    if not empty.empty:
        z = pd.DataFrame({"locus_id": empty.index})
        for col in ("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"):
            z[col] = pd.NA
        res = pd.concat([res, z], ignore_index=True)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(out_tsv, sep="\t", index=False)
    if len(res) != len(counts):
        log.warning("de.tsv rows %s != locus_counts %s", len(res), len(counts))
    return res


DE_STATUS_WALD = "wald"
DE_STATUS_LOW_COUNT = "low_count"


def de_status_series(res: pd.DataFrame, cfg: TxNovaConfig) -> pd.Series:
    """Label each DE row. Empty string = do not keep.

    `wald`: padj and LFC both meet the configured gates.
    `low_count`: DESeq2 independent filtering left padj NA but Wald still
    ran (pvalue present) and LFC is treat-up. That is "not evaluated",
    not a reject and not a significant padj.
    Cook / all-zero / not fitted (pvalue also NA) stay empty.
    Evaluated non-hits (padj present but too large, or LFC too small) stay empty.
    """
    if res.empty:
        return pd.Series(dtype=object)
    if "pvalue" not in res.columns:
        raise TxNovaError("de.tsv missing column pvalue; delete stamps/de.json and rerun")
    padj = pd.to_numeric(res["padj"], errors="coerce")
    lfc = pd.to_numeric(res["log2FoldChange"], errors="coerce")
    pvalue = pd.to_numeric(res["pvalue"], errors="coerce")
    up = lfc.notna() & (lfc >= cfg.de.min_log2fc)
    wald = padj.notna() & (padj < cfg.de.padj) & up
    low_count = padj.isna() & pvalue.notna() & up
    status = pd.Series("", index=res.index, dtype=object)
    status = status.mask(low_count, DE_STATUS_LOW_COUNT)
    status = status.mask(wald, DE_STATUS_WALD)
    return status


def de_pass_mask(res: pd.DataFrame, cfg: TxNovaConfig) -> pd.Series:
    status = de_status_series(res, cfg)
    return status.isin([DE_STATUS_WALD, DE_STATUS_LOW_COUNT])
