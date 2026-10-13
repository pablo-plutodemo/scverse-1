# Human white adipose tissue snRNA-seq: depot comparison of adipocytes

Reproducible analysis of the CELLxGENE dataset *Single-Nucleus Analysis of Human
White Adipose Tissue Reveals Adipocyte Subsets with Distinct Metabolic Profiles*
(collection `6b701826-37bb-4356-9792-ff41fc4c3161`, 65,668 nuclei, 19 donors,
subcutaneous and visceral abdominal depots).

## What the pipeline does

`notebooks/wat_depot_analysis.ipynb` runs end to end:

1. **Download / validate** the H5AD from CELLxGENE (skipped if the file already
   exists; schema, counts and normalisation are checked instead).
2. **Overview UMAPs** coloured by *ADIPOQ* and by cell type.
3. **Identify the "unknown" cluster** with a donor-stratified t-test against every
   other cell type (Stouffer-combined). It is a spermatid transcript signature
   (*PRM1*, *TNP1*, *ODF2*, *TUBA3C*, …) present mostly in female donors, i.e. a
   technical contaminant, and is excluded downstream.
4. **Build the adipocyte-only object** (raw UMI counts, depot/donor/library fields).
5. **Visceral vs subcutaneous adipocyte DE** — pseudobulk per library with
   pyDESeq2 (`~ depot`, 12 vs 10 libraries) as the primary analysis, and a
   cell-level donor-stratified t-test (restricted to the 3 donors with both
   depots) alongside for comparison.
6. **Volcano plots** and a full per-gene table.

## Running

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
export WAT_DATA_DIR=~/adipose_cxg/data      # where the ~1 GB h5ad lives / will be saved
jupyter nbconvert --to notebook --execute --inplace notebooks/wat_depot_analysis.ipynb
```

Outputs are written to `figures/` and `results/`. Peak memory is ~12 GB.

## Layout

```
notebooks/wat_depot_analysis.ipynb   the pipeline
src/wat_pipeline.py                  helper functions used by the notebook
figures/                             PNG outputs (committed)
results/                             CSV tables (git-ignored except small summaries)
requirements.txt                     frozen dependencies
```

## Notes

- The dataset's `X` is log-normalised; raw integer counts are in `adata.raw.X`.
  `scanpy.tl.rank_genes_groups` defaults to `adata.raw` here, so `use_raw=False`
  is passed explicitly.
- Depot is nested in donor (only 3 of 19 donors contribute both depots), so a
  cell-level "donor-adjusted" test can only use those paired donors; the
  pseudobulk result is the one to cite.
- pyDESeq2 is run with the joblib threading backend (see `run_deseq2`) because
  the default loky backend needs POSIX semaphores that sandboxed macOS forbids.
