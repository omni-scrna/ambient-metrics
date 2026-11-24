#!/usr/bin/env python3
"""Merged vs. unmerged input, for one ambient removal method (stage MERGE_M).

ambient-plan runs every method on the raw counts (PREP `merge: none`) and, where
it applies, on merged counts (`orthologs` for mixed-species pools, `gems` for
Flex pools). MERGE_M is a gather with `split_by: AMBIENT`: one node per method
and parameter set, whose members are that method's outputs on every PREP
variant and pool.

For each pool, the unmerged output is folded into the merged id space (the same
sums PREP applies, taken from the names: cells by GEM barcode, genes by ortholog
pair) and compared with the merged output. Does removing ambient RNA and then
merging give what merging first gives?

Members come from the manifest (`--from-files lineage.json`, omnibenchmark
design 010 §3.5.1): the PREP variant is the `merge` parameter of a member's
PREP lineage entry, the pool is its `value`. Writes one JSON with a record per
pair; a dataset with no merged variant (HTO) gives `n_pairs: 0`.
"""

import argparse
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))
from common import cli  # noqa: E402  (src/common, copied from boilerplate)
from run import N_ALIGN_CELLS, _cor, _f, _ratio, _rowcor, strip_output_args  # noqa: E402

OUTPUT_ID = "merge_metrics_json"
FROM_ID = "filtered_h5ad"
VARIANT_STAGE = "PREP"
UNMERGED = "none"


def gem_of(cell_id):
    """'{gem_barcode}_{probe_id}' -> '{gem_barcode}' (PREP `merge: gems`)."""
    return cell_id.rsplit("_", 1)[0]


def bare_gene(name):
    """'GRCh38_ENSG…' -> 'ENSG…': the id without its genome prefix."""
    return name.split("_", 1)[1] if "_" in name else name


def _indicator(source, target_of, targets):
    """Sparse len(source) x len(targets) 0/1 matrix; sources without a target get no column."""
    col = {t: j for j, t in enumerate(targets)}
    rows, cols = [], []
    for i, s in enumerate(source):
        j = col.get(target_of(s))
        if j is not None:
            rows.append(i)
            cols.append(j)
    return sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(source), len(targets)))


def fold(unmerged, merged):
    """`unmerged` counts summed into the cells x genes of `merged`.

    Names decide the sums, as PREP's merges do: a cell missing from `merged`
    maps to its GEM barcode, a gene to the ortholog pair that names it
    ('{human_id}_{mouse_id}'). Unmapped cells and genes are dropped.
    """
    cells_m, genes_m = set(merged.obs_names), set(merged.var_names)
    cell_target = (lambda c: c) if set(unmerged.obs_names) <= cells_m else gem_of
    if set(unmerged.var_names) <= genes_m:
        gene_target = lambda g: g  # noqa: E731
    else:
        pair_of = {}
        for pair in merged.var_names:
            for part in pair.split("_"):
                pair_of[part] = pair
        gene_target = lambda g: pair_of.get(bare_gene(g))  # noqa: E731
    C = _indicator(unmerged.obs_names, cell_target, list(merged.obs_names))
    G = _indicator(unmerged.var_names, gene_target, list(merged.var_names))
    X = sp.csr_matrix(unmerged.X).astype(np.float64)
    return (C.T @ X @ G).tocsr(), int(C.sum()), int(G.sum())


def compare_pair(unmerged, merged, seed=0):
    """Metrics for one method on one pool: unmerged output folded vs. merged output."""
    F, n_cells_mapped, n_genes_mapped = fold(unmerged, merged)
    M = sp.csr_matrix(merged.X).astype(np.float64)
    # cells the unmerged run also has (a GEM-merged run may call the same GEMs only)
    rows = np.where(np.asarray(F.sum(axis=1)).ravel() + np.asarray(M.sum(axis=1)).ravel() > 0)[0]
    F, M = F[rows], M[rows]
    tot_f, tot_m = float(F.sum()), float(M.sum())
    m = {
        "n_cells_unmerged": int(unmerged.n_obs),
        "n_cells_merged": int(merged.n_obs),
        "n_cells_unmerged_mapped": n_cells_mapped,
        "n_genes_unmerged": int(unmerged.n_vars),
        "n_genes_merged": int(merged.n_vars),
        "n_genes_unmerged_mapped": n_genes_mapped,
        "total_counts_unmerged_folded": tot_f,
        "total_counts_merged": tot_m,
        "rel_total_diff": _f((tot_m - tot_f) / tot_f) if tot_f else None,
        "rel_l1_diff": _ratio(float(abs(F - M).sum()), tot_f),
    }
    n = F.shape[0]
    if n >= 2:
        idx = np.random.default_rng(seed).choice(n, size=min(N_ALIGN_CELLS, n), replace=False)
        m["cell_corr_median"] = _f(np.nanmedian(_rowcor(F[idx], M[idx])))
    else:
        m["cell_corr_median"] = None
    m["gene_corr"] = _f(_cor(np.asarray(F.sum(axis=0)).ravel(), np.asarray(M.sum(axis=0)).ravel()))
    return m


def variant_of(member):
    """The PREP variant a member ran on: the `merge` parameter, else the module id."""
    for entry in member["lineage"]:
        if entry["stage"] == VARIANT_STAGE:
            return str((entry.get("parameters") or {}).get("merge", entry["module"]))
    raise ValueError(f"member {member['node_id']} has no {VARIANT_STAGE} entry in its lineage")


def pairs_of(manifest):
    """[(pool, variant, unmerged member, merged member)], one per merged member."""
    by_pool = {}
    for member in manifest["members"]:
        by_pool.setdefault(member.get("value"), {})[variant_of(member)] = member
    out = []
    for pool, variants in sorted(by_pool.items(), key=lambda kv: str(kv[0])):
        unmerged = variants.get(UNMERGED)
        if unmerged is None:
            continue
        out += [(pool, v, unmerged, m) for v, m in sorted(variants.items()) if v != UNMERGED]
    return out


def compare_manifest(manifest, seed=0):
    method = manifest["members"][0]["lineage"][-1] if manifest["members"] else {}
    records = []
    for pool, variant, a, b in pairs_of(manifest):
        rec = {"pool": pool, "merge": variant}
        rec.update(compare_pair(ad.read_h5ad(a["files"][FROM_ID]), ad.read_h5ad(b["files"][FROM_ID]), seed))
        records.append(rec)
    return {
        "method": method.get("module"),
        "parameters": method.get("parameters") or {},
        "n_pairs": len(records),
        "pairs": records,
    }


def main():
    rest, outputs = strip_output_args(sys.argv[1:])
    p = argparse.ArgumentParser(description="Merged vs. unmerged input, per method", allow_abbrev=False)
    cli.add_base_args(p)               # --output_dir, --name
    cli.add_stage_args(p, "MERGE_M")   # --from-files
    p.add_argument("--seed", type=int, default=0, help="Seed for the cell sample of the correlation")
    args = p.parse_args(rest)
    print(f"Full command: {' '.join(sys.argv)}")

    with open(args.from_files) as fh:
        manifest = json.load(fh)
    result = compare_manifest(manifest, seed=args.seed)
    path = Path(outputs.get(OUTPUT_ID) or Path(args.output_dir) / f"{args.name}_merge_metrics.json")
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(result, fh, indent=2)
    print(f"  wrote: {path} ({result['n_pairs']} pair(s))")


if __name__ == "__main__":
    main()
