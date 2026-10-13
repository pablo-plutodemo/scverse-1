"""Helper functions for the WAT snRNA-seq depot analysis notebook."""
from __future__ import annotations

import os
import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import requests
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

COLLECTION_ID = "6b701826-37bb-4356-9792-ff41fc4c3161"
DATASET_ID = "9d8e5dca-03a3-457d-b7fb-844c75735c83"
CXG_API = "https://api.cellxgene.cziscience.com/curation/v1"


def data_dir() -> Path:
    d = Path(os.environ.get("WAT_DATA_DIR", "~/adipose_cxg/data")).expanduser()
    d.mkdir(parents=True, exist_ok=True)
    return d


def download_or_validate(path: Path) -> ad.AnnData:
    """Download the H5AD from CELLxGENE unless it exists; always validate and return a backed handle."""
    if not path.exists():
        meta = requests.get(f"{CXG_API}/collections/{COLLECTION_ID}/datasets/{DATASET_ID}", timeout=60).json()
        url = [a["url"] for a in meta["assets"] if a["filetype"] == "H5AD"][0]
        tmp = path.with_suffix(".h5ad.part")
        with requests.get(url, stream=True, timeout=600) as r:
            r.raise_for_status()
            with open(tmp, "wb") as f:
                for chunk in r.iter_content(1 << 22):
                    f.write(chunk)
        tmp.rename(path)
    a = ad.read_h5ad(path, backed="r")
    assert a.raw is not None, "expected raw counts in adata.raw"
    assert "cell_type" in a.obs and "donor_id" in a.obs and "X_umap" in a.obsm
    probe = a[:1000].to_memory()
    assert np.allclose(probe.raw.X.data, np.round(probe.raw.X.data)), "raw.X should be integer UMI counts"
    assert probe.X.data.max() < 30, "X should be log-normalised"
    return a


def gene_symbols(a: ad.AnnData) -> None:
    """Switch var_names from Ensembl IDs to feature_name symbols (in place)."""
    a.var_names = a.var["feature_name"].astype(str).values
    a.var_names_make_unique()


def stratified_ttest(adata: ad.AnnData, group_key: str, group: str, reference: str,
                     stratum_key: str, min_cells: int = 20) -> pd.DataFrame:
    """Per-stratum scanpy t-test (group vs reference) combined across strata with Stouffer's method.

    Weights are sqrt(n cells in stratum). Returns z, p, BH padj, log2FC of expm1 means,
    detection fractions and the fraction of strata whose sign agrees with the combined z.
    """
    ts, ws, used = [], [], []
    for s in adata.obs[stratum_key].unique():
        m = (adata.obs[stratum_key] == s) & adata.obs[group_key].isin([group, reference])
        sub = adata[m]
        n_g = (sub.obs[group_key] == group).sum(); n_r = (sub.obs[group_key] == reference).sum()
        if n_g < min_cells or n_r < min_cells:
            continue
        sub = sub.copy(); sub.obs["_g"] = sub.obs[group_key].astype(str)
        sc.tl.rank_genes_groups(sub, "_g", groups=[group], reference=reference, method="t-test",
                                n_genes=sub.n_vars, use_raw=False)
        df = sc.get.rank_genes_groups_df(sub, group).set_index("names")
        ts.append(df["scores"].reindex(adata.var_names)); ws.append(np.sqrt(sub.n_obs)); used.append(s)
    if not ts:
        return pd.DataFrame()
    T = pd.concat(ts, axis=1).fillna(0); w = np.array(ws)
    z = (T * w).sum(1) / np.sqrt((w ** 2).sum())
    p = 2 * norm.sf(np.abs(z))
    m = adata.obs[stratum_key].isin(used)
    G = adata[m & (adata.obs[group_key] == group)].X; R = adata[m & (adata.obs[group_key] == reference)].X
    mu_g = np.asarray(np.expm1(G).mean(0)).ravel(); mu_r = np.asarray(np.expm1(R).mean(0)).ravel()
    return pd.DataFrame({
        "stouffer_z": z.values, "pval": p, "padj": multipletests(p, method="fdr_bh")[1],
        "log2FC": np.log2((mu_g + 1e-2) / (mu_r + 1e-2)),
        "pct_group": np.asarray((G > 0).mean(0)).ravel(),
        "pct_reference": np.asarray((R > 0).mean(0)).ravel(),
        "sign_consistency": (np.sign(T).T == np.sign(z)).mean(0).values,
        "n_strata": len(used), "n_group": G.shape[0], "n_reference": R.shape[0],
    }, index=adata.var_names)


def build_adipocyte_object(a: ad.AnnData, out: Path) -> ad.AnnData:
    """Adipocyte-only AnnData with raw UMI counts in X; load from `out` if it exists."""
    if out.exists():
        adi = ad.read_h5ad(out)
    else:
        m = (a.obs["cell_type"] == "adipocyte").values
        adi = ad.AnnData(X=a.raw.X[m].astype(np.int32), obs=a.obs.loc[m].copy(), var=a.raw.var.copy())
        adi.var["feature_name"] = a.var["feature_name"].values
        adi.var_names = adi.var["feature_name"].astype(str).values; adi.var_names_make_unique()
        adi.obs["depot"] = adi.obs["origin"].map({"SQ": "subcutaneous", "GC": "visceral", "OM": "visceral"}).astype("category")
        adi.obs["donor"] = adi.obs["donor_id"].astype(str).str.replace("donor-", "").astype("category")
        adi.obs["library"] = adi.obs["library_id"].astype(str).astype("category")
        adi.obsm["X_umap"] = a.obsm["X_umap"][m]
        adi.uns["source"] = {"file": "wat_snrna_original.h5ad", "cellxgene_dataset_id": DATASET_ID,
                             "selection": "cell_type == 'adipocyte'; 'unknown' (PRM1+ spermatid contaminant) excluded",
                             "X": "raw integer UMI counts from original adata.raw.X",
                             "created": datetime.date.today().isoformat()}
        adi.write_h5ad(out)
    assert adi.X.dtype.kind == "i" and adi.X.max() > 50, "expected raw integer counts"
    assert {"depot", "donor", "library"} <= set(adi.obs.columns)
    return adi


def pseudobulk(adi: ad.AnnData, by: str = "library") -> tuple[pd.DataFrame, pd.DataFrame]:
    """Sum raw counts per `by` level; return (counts [samples x genes], sample metadata)."""
    levels = adi.obs[by].cat.categories
    counts = pd.DataFrame(np.vstack([np.asarray(adi[adi.obs[by] == l].X.sum(0)).ravel() for l in levels]),
                          index=levels, columns=adi.var_names)
    meta = adi.obs.groupby(by, observed=True).agg(depot=("depot", "first"), donor=("donor", "first"),
                                                  T2D=("T2D", "first"), sex=("sex", "first"),
                                                  n_cells=("depot", "size")).loc[levels]
    keep = ((counts > 0).sum(0) >= 3) & (counts.sum(0) >= 10)
    return counts.loc[:, keep], meta


def run_deseq2(counts: pd.DataFrame, meta: pd.DataFrame, contrast=("depot", "visceral", "subcutaneous"),
               n_cpus: int = 4) -> pd.DataFrame:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    from pydeseq2.default_inference import DefaultInference
    # pyDESeq2 wraps joblib with inner_max_num_threads=1, which only the loky backend accepts; loky in turn
    # needs POSIX semaphores that sandboxed macOS forbids. Use the threading backend and tell joblib it is fine.
    from joblib._parallel_backends import ThreadingBackend
    ThreadingBackend.supports_inner_max_num_threads = True
    inf = DefaultInference(n_cpus=n_cpus, backend="threading", joblib_verbosity=0)
    dds = DeseqDataSet(counts=counts, metadata=meta, design=f"~{contrast[0]}", quiet=True, inference=inf)
    dds.deseq2()
    st = DeseqStats(dds, contrast=list(contrast), quiet=True, inference=inf)
    st.summary()
    return st.results_df.rename(columns={"log2FoldChange": "pb_log2FC", "lfcSE": "pb_lfcSE", "stat": "pb_stat",
                                         "pvalue": "pb_pval", "padj": "pb_padj", "baseMean": "pb_baseMean"})
