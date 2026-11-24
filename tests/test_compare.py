import json
import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from compare import compare_manifest, compare_pair, fold, pairs_of  # noqa: E402


def adata(X, cells, genes):
    return ad.AnnData(
        X=csr_matrix(np.asarray(X, dtype=np.int64)),
        obs=pd.DataFrame(index=cells),
        var=pd.DataFrame(index=genes),
    )


# Mixed species: two 1:1 ortholog pairs, plus a human gene with no partner.
HUMAN_MOUSE = ["GRCh38_ENSG1", "GRCh38_ENSG2", "GRCh38_ENSG9", "GRCm39_ENSMUSG1", "GRCm39_ENSMUSG2"]
PAIRS = ["ENSG1_ENSMUSG1", "ENSG2_ENSMUSG2"]


def test_fold_sums_ortholog_pairs():
    unmerged = adata([[1, 2, 50, 10, 20], [3, 4, 60, 30, 40]], ["c1", "c2"], HUMAN_MOUSE)
    merged = adata([[11, 22], [33, 44]], ["c1", "c2"], PAIRS)
    F, n_cells, n_genes = fold(unmerged, merged)
    assert F.toarray().tolist() == [[11, 22], [33, 44]]
    assert (n_cells, n_genes) == (2, 4)  # ENSG9 has no pair
    m = compare_pair(unmerged, merged)
    assert m["rel_l1_diff"] == 0 and m["rel_total_diff"] == 0 and m["gene_corr"] > 0.99


def test_fold_sums_gems():
    # G1 was read under two probes; the merged run has one cell per GEM, reordered.
    unmerged = adata([[1, 2], [10, 20], [5, 5]], ["G1_BC001", "G1_BC002", "G2_BC001"], ["g1", "g2"])
    merged = adata([[5, 5], [11, 22]], ["G2", "G1"], ["g1", "g2"])
    F, n_cells, _ = fold(unmerged, merged)
    assert F.toarray().tolist() == [[5, 5], [11, 22]]
    assert n_cells == 3


def test_difference_is_measured():
    unmerged = adata([[10, 0], [0, 10]], ["c1", "c2"], ["g1", "g2"])
    merged = adata([[8, 0], [0, 10]], ["c1", "c2"], ["g1", "g2"])
    m = compare_pair(unmerged, merged)
    assert m["rel_total_diff"] == -0.1 and m["rel_l1_diff"] == 0.1


def member(tmp_path, pool, merge, X, cells, genes, method="decontx"):
    path = tmp_path / f"{pool}_{merge}_{method}.h5ad"
    adata(X, cells, genes).write_h5ad(path)
    return {
        "node_id": f"n-{merge}-{method}",
        "value": pool,
        "files": {"filtered_h5ad": str(path)},
        "lineage": [
            {"stage": "DATA_RAW", "module": "d_flex", "parameters": {}},
            {"stage": "PREP", "module": f"prep_{merge}", "parameters": {"merge": merge}},
            {"stage": "AMBIENT", "module": method, "parameters": {"use_background": "true"}},
        ],
    }


def manifest(tmp_path):
    cells, genes = ["G1_BC001", "G1_BC002"], ["g1"]
    return {"members": [
        member(tmp_path, "p1", "none", [[1], [2]], cells, genes),
        member(tmp_path, "p2", "none", [[4], [4]], cells, genes),
        member(tmp_path, "p1", "gems", [[3]], ["G1"], genes),
        member(tmp_path, "p2", "gems", [[6]], ["G1"], genes),
    ]}


def test_pairs_by_pool_from_the_manifest(tmp_path):
    man = manifest(tmp_path)
    assert [(pool, v) for pool, v, _a, _b in pairs_of(man)] == [("p1", "gems"), ("p2", "gems")]
    res = compare_manifest(man)
    assert res["method"] == "decontx" and res["parameters"] == {"use_background": "true"}
    assert [(r["pool"], r["rel_total_diff"]) for r in res["pairs"]] == [("p1", 0.0), ("p2", -0.25)]


def test_no_merged_variant_gives_no_pairs(tmp_path):
    man = {"members": [member(tmp_path, "p1", "none", [[1]], ["c1"], ["g1"])]}
    assert compare_manifest(man) == {
        "method": "decontx", "parameters": {"use_background": "true"}, "n_pairs": 0, "pairs": []}


def test_cli(tmp_path):
    (tmp_path / "lineage.json").write_text(json.dumps(manifest(tmp_path)))
    out = tmp_path / "out" / "m.json"
    subprocess.run(
        [sys.executable, str(ROOT / "compare.py"), "--output_dir", str(tmp_path), "--name", "cmp",
         "--from-files", str(tmp_path / "lineage.json"), "--output", f"merge_metrics_json={out}"],
        check=True, capture_output=True)
    assert json.loads(out.read_text())["n_pairs"] == 2
