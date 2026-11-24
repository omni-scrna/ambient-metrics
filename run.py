#!/usr/bin/env python3
"""Ambient RNA removal metrics (stage AMBIENT_M of ambient-plan).

Compares one ambient-removal output (`filtered_h5ad`) of one pool with the
counts the method was given for the same cells (`prepared_h5ad`, the PREP
output: raw, or merged over orthologs / GEMs), and writes one flat JSON of
scalars. "Raw" below means that matrix.

1. Count inflation (every dataset). Ambient removal should only subtract
   counts. Does a method *add* counts compared with the raw counts?
2. Barnyard (mixed-species pools only, e.g. HGMM). In a human cell, counts on
   mouse genes can only come from ambient RNA (and vice versa). The species of
   each cell comes from the DEMUX stage (run on the raw counts), the species of
   each gene from the raw matrix (`var['species']`). Metrics:
   - cross-species fraction of the counts before and after removal;
   - fraction of cross-species counts removed (ideal: close to 1);
   - fraction of own-species counts removed (should stay small: own-species
     ambient RNA exists too, so it is not 0 for a perfect method).
   Datasets without both species give None for these keys, and so does an
   ortholog-merged matrix: its features belong to neither species.

Everything is aligned by barcode and gene NAME, never by position: methods
may reorder cells and genes. A name-independent check (correlation of matched
vs. shifted cells) catches misaligned names. `no_cleaning` is the control:
0 counts above raw, 0 removed, alignment ~1.
"""

import argparse
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))
from common import cli  # noqa: E402  (src/common, copied from boilerplate)

OUTPUT_ID = "ambient_metrics_json"
N_ALIGN_CELLS = 500
SPECIES = ("human", "mouse")


# ---------------------------------------------------------------- CLI

def strip_output_args(argv):
    """api >= 0.8 passes `--output <id>=<path>` once per output. The boilerplate
    CLI doesn't declare it, so take those pairs out before parsing."""
    rest, outputs, i = [], {}, 0
    while i < len(argv):
        if argv[i] == "--output" and i + 1 < len(argv):
            key, _, path = argv[i + 1].partition("=")
            outputs[key] = path
            i += 2
        elif argv[i].startswith("--output="):
            key, _, path = argv[i].split("=", 1)[1].partition("=")
            outputs[key] = path
            i += 1
        else:
            rest.append(argv[i])
            i += 1
    return rest, outputs


def parse_args(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    rest, outputs = strip_output_args(argv)
    p = argparse.ArgumentParser(description="Ambient RNA removal metrics", allow_abbrev=False)
    cli.add_base_args(p)                  # --output_dir, --name
    cli.add_stage_args(p, "AMBIENT_M")    # --prepared_h5ad, --filtered_h5ad, --demultiplexed_tsv
    p.add_argument("--seed", type=int, default=0, help="Seed for the alignment-check cell sample")
    args = p.parse_args(rest)
    args.outputs = outputs
    return args


def output_path(args):
    if OUTPUT_ID in args.outputs:
        return Path(args.outputs[OUTPUT_ID])
    return Path(args.output_dir) / f"{args.name}_ambient_metrics.json"


# ---------------------------------------------------------------- helpers

def _f(x):
    x = float(x)
    return None if np.isnan(x) else x


def _ratio(a, b):
    return _f(a / b) if b else None


def _rowcor(a, b):
    """Pearson correlation of log1p counts, row by row (dense; few rows)."""
    a = np.log1p(np.asarray(a.todense(), dtype=float))
    b = np.log1p(np.asarray(b.todense(), dtype=float))
    a -= a.mean(axis=1, keepdims=True)
    b -= b.mean(axis=1, keepdims=True)
    den = np.sqrt((a**2).sum(axis=1) * (b**2).sum(axis=1))
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, (a * b).sum(axis=1) / den, np.nan)


def _cor(x, y):
    x = np.log1p(np.asarray(x, dtype=float))
    y = np.log1p(np.asarray(y, dtype=float))
    if x.std() == 0 or y.std() == 0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def _rowsums(M, cols):
    return np.asarray(M[:, cols].sum(axis=1)).ravel().astype(float)


# ---------------------------------------------------------------- metrics

def align(raw, out):
    """Validity checks; return (metrics, R, Xf, out) aligned by name.
    R: raw counts (int64), Xf: output counts (float, unrounded), same cells x genes."""
    m = {"n_cells": int(out.n_obs), "n_genes": int(out.n_vars)}
    if not raw.obs_names.is_unique or not raw.var_names.is_unique:
        raise ValueError("raw barcodes/genes are not unique; cannot align by name")
    if not out.obs_names.is_unique or not out.var_names.is_unique:
        raise ValueError("output barcodes/genes are not unique; cannot align by name")
    cells_in_raw = out.obs_names.isin(raw.obs_names)
    genes_in_raw = out.var_names.isin(raw.var_names)
    m["n_cells_not_in_raw"] = int((~cells_in_raw).sum())
    m["n_genes_not_in_raw"] = int((~genes_in_raw).sum())
    out = out[cells_in_raw, genes_in_raw]
    Xf = sp.csr_matrix(out.X).astype(np.float64)
    m["is_integer"] = bool(np.all(np.mod(Xf.data, 1) == 0))
    m["n_negative"] = int((Xf.data < 0).sum())
    R = sp.csr_matrix(raw[out.obs_names, out.var_names].X).astype(np.int64)
    return m, R, Xf, out


def inflation_metrics(R, Xf, seed=0):
    m = {}
    X = Xf.copy()
    X.data = np.round(X.data)
    X = X.astype(np.int64)

    tot_raw, tot_out = int(R.sum()), float(Xf.sum())
    m["total_counts_raw"] = tot_raw
    m["total_counts_out"] = tot_out
    m["frac_counts_removed"] = _f(1 - tot_out / tot_raw) if tot_raw else None

    D = (X - R).tocsr()
    above = D.multiply(D > 0).tocsr()
    above.eliminate_zeros()
    excess = int(above.sum())
    m["n_entries_above_raw"] = int(above.nnz)
    m["frac_entries_above_raw"] = _ratio(above.nnz, R.nnz)
    m["excess_counts"] = excess
    m["frac_excess_of_out"] = _ratio(excess, tot_out)
    # 0 in raw, > 0 in output; avoid `R == 0` (that matrix is almost dense)
    X_pos, R_pos = (X > 0).astype(np.int8), (R > 0).astype(np.int8)
    m["n_new_nonzeros"] = int(X_pos.nnz - X_pos.multiply(R_pos).nnz)

    cell_excess = np.asarray(above.sum(axis=1)).ravel()
    cell_out = np.asarray(X.sum(axis=1)).ravel()
    with np.errstate(invalid="ignore", divide="ignore"):
        cell_frac = np.where(cell_out > 0, cell_excess / cell_out, 0.0)
    has = len(cell_frac) > 0
    m["n_cells_with_excess"] = int((cell_excess > 0).sum())
    m["frac_cells_with_excess"] = _f((cell_excess > 0).mean()) if has else None
    m["cell_excess_frac_median"] = _f(np.median(cell_frac)) if has else None
    m["cell_excess_frac_p90"] = _f(np.quantile(cell_frac, 0.9)) if has else None
    m["cell_excess_frac_max"] = _f(cell_frac.max()) if has else None
    m["n_genes_with_excess"] = int((np.asarray(above.sum(axis=0)).ravel() > 0).sum())

    n = X.shape[0]
    if n >= 2:
        rng = np.random.default_rng(seed)
        idx = rng.choice(n, size=min(N_ALIGN_CELLS, n), replace=False)
        m["align_cell_corr_matched"] = _f(np.nanmedian(_rowcor(R[idx], X[idx])))
        m["align_cell_corr_shifted"] = _f(np.nanmedian(_rowcor(R[np.roll(idx, 1)], X[idx])))
    else:
        m["align_cell_corr_matched"] = m["align_cell_corr_shifted"] = None
    g_raw = np.asarray(R.sum(axis=0)).ravel()
    g_out = np.asarray(X.sum(axis=0)).ravel()
    m["align_gene_corr_matched"] = _f(_cor(g_raw, g_out))
    m["align_gene_corr_shifted"] = _f(_cor(np.roll(g_raw, 1), g_out))
    return m


BARNYARD_KEYS = (
    ["barnyard", "n_singlets", "n_cells_no_label", "n_multiplets_skipped"]
    + [f"{s}_{k}" for s in ("all",) + SPECIES for k in (
        "n_cells", "xspecies_frac_raw", "xspecies_frac_out",
        "xspecies_removed", "own_removed",
        "cell_xspecies_frac_raw_median", "cell_xspecies_frac_out_median")]
)


def barnyard_metrics(R, Xf, cells, gene_species, demux):
    """Cross-species metrics. `cells`: obs names of R/Xf; `gene_species`: array of
    'human'/'mouse'/other per column; `demux`: DataFrame with cell_id, species,
    multiplet (DEMUX output)."""
    m = dict.fromkeys(BARNYARD_KEYS)
    m["barnyard"] = False
    gene_species = np.asarray(gene_species, dtype=object)
    if demux is None or "species" not in demux.columns:
        return m
    if not all((gene_species == s).any() for s in SPECIES):
        return m

    d = demux.drop_duplicates("cell_id").set_index("cell_id").reindex(cells)
    multiplet = d["multiplet"].astype(str).str.upper().eq("TRUE").to_numpy() \
        if "multiplet" in d.columns else np.zeros(len(cells), dtype=bool)
    cell_species = d["species"].to_numpy(dtype=object)
    labelled = pd.notna(d["species"]).to_numpy()
    m["barnyard"] = True
    m["n_cells_no_label"] = int((~labelled).sum())
    m["n_multiplets_skipped"] = int((labelled & multiplet).sum())

    sums = {}
    for s in SPECIES:
        other = [o for o in SPECIES if o != s][0]
        rows = np.where(labelled & ~multiplet & (cell_species == s))[0]
        own_cols = np.where(gene_species == s)[0]
        oth_cols = np.where(gene_species == other)[0]
        Rs, Xs = R[rows], Xf[rows]
        sums[s] = dict(
            n=len(rows),
            raw_own=_rowsums(Rs, own_cols), raw_oth=_rowsums(Rs, oth_cols),
            out_own=_rowsums(Xs, own_cols), out_oth=_rowsums(Xs, oth_cols),
        )
    m["n_singlets"] = int(sum(v["n"] for v in sums.values()))

    groups = {s: sums[s] for s in SPECIES}
    groups["all"] = {k: (sum(sums[s][k] for s in SPECIES) if k == "n"
                         else np.concatenate([sums[s][k] for s in SPECIES]))
                     for k in sums[SPECIES[0]]}
    for g, v in groups.items():
        ro, rx, oo, ox = (v[k].sum() for k in ("raw_own", "raw_oth", "out_own", "out_oth"))
        m[f"{g}_n_cells"] = int(v["n"])
        m[f"{g}_xspecies_frac_raw"] = _ratio(rx, ro + rx)
        m[f"{g}_xspecies_frac_out"] = _ratio(ox, oo + ox)
        m[f"{g}_xspecies_removed"] = _f(1 - ox / rx) if rx else None
        m[f"{g}_own_removed"] = _f(1 - oo / ro) if ro else None
        with np.errstate(invalid="ignore", divide="ignore"):
            fr = v["raw_oth"] / (v["raw_own"] + v["raw_oth"])
            fo = v["out_oth"] / (v["out_own"] + v["out_oth"])
        m[f"{g}_cell_xspecies_frac_raw_median"] = _f(np.nanmedian(fr)) if np.isfinite(fr).any() else None
        m[f"{g}_cell_xspecies_frac_out_median"] = _f(np.nanmedian(fo)) if np.isfinite(fo).any() else None
    return m


def compute_metrics(raw, out, demux=None, seed=0):
    """All metrics for one method output `out` vs. the raw droplets `raw` (AnnData)."""
    m, R, Xf, out = align(raw, out)
    m.update(inflation_metrics(R, Xf, seed=seed))
    if "species" in raw.var.columns:
        gene_species = raw.var["species"].reindex(out.var_names).to_numpy(dtype=object)
    else:
        gene_species = np.full(out.n_vars, None, dtype=object)
    m.update(barnyard_metrics(R, Xf, out.obs_names, gene_species, demux))
    return m


# ---------------------------------------------------------------- main

def main():
    args = parse_args()
    print(f"Full command: {' '.join(sys.argv)}")
    for k, v in vars(args).items():
        print(f"  {k}: {v}")

    # backed="r": a raw pool has ~1.5 M barcodes; load only the output's cells.
    raw = ad.read_h5ad(args.prepared_h5ad, backed="r")
    out = ad.read_h5ad(args.filtered_h5ad)
    keep = out.obs_names[out.obs_names.isin(raw.obs_names)]
    raw = raw[keep].to_memory()
    demux = pd.read_csv(args.demultiplexed_tsv, sep="\t", dtype={"cell_id": str})
    print(f"  output: {out.n_obs} cells x {out.n_vars} genes; raw subset: {raw.n_obs} cells; "
          f"demux: {len(demux)} rows")

    metrics = compute_metrics(raw, out, demux=demux, seed=args.seed)
    path = output_path(args)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(metrics, fh, indent=2)
    print(f"  wrote: {path}")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
