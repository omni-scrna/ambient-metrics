import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from run import compute_metrics, output_path, parse_args  # noqa: E402


def _adata(X, cells, genes, species=None):
    var = pd.DataFrame(index=genes)
    if species is not None:
        var["species"] = species
    return ad.AnnData(X=csr_matrix(np.asarray(X)), obs=pd.DataFrame(index=cells), var=var)


@pytest.fixture
def raw():
    rng = np.random.default_rng(1)
    X = rng.poisson(3, size=(40, 30))
    return _adata(X, [f"c{i}" for i in range(40)], [f"g{j}" for j in range(30)])


# ---------------------------------------------------------------- count inflation

def test_passthrough_is_clean(raw):
    out = raw[:20].copy()  # "no_cleaning": raw counts of the called cells
    m = compute_metrics(raw, out)
    assert m["n_cells_not_in_raw"] == 0 and m["n_genes_not_in_raw"] == 0
    assert m["n_entries_above_raw"] == 0 and m["excess_counts"] == 0
    assert m["n_new_nonzeros"] == 0
    assert m["frac_counts_removed"] == pytest.approx(0)
    assert m["align_cell_corr_matched"] == pytest.approx(1)
    assert m["align_gene_corr_matched"] == pytest.approx(1)
    assert m["barnyard"] is False and m["all_xspecies_removed"] is None


def test_reordered_output_is_aligned_by_name(raw):
    out = raw[:20].copy()
    out = out[out.obs_names[::-1], out.var_names[::-1]].copy()
    m = compute_metrics(raw, out)
    assert m["n_entries_above_raw"] == 0 and m["excess_counts"] == 0
    assert m["align_cell_corr_matched"] == pytest.approx(1)


def test_index_mismatch_would_show_as_inflation(raw):
    sub = raw[:20].copy()
    shuffled = _adata(sub.X.toarray()[::-1], list(sub.obs_names), list(sub.var_names))
    m = compute_metrics(raw, shuffled)
    assert m["n_entries_above_raw"] > 0
    assert m["align_cell_corr_matched"] < 0.5


def test_counts_added(raw):
    sub = raw[:20].copy()
    X = sub.X.toarray()
    X[0, 0] += 5
    X[1, :] = 0
    out = _adata(X, list(sub.obs_names), list(sub.var_names))
    m = compute_metrics(raw, out)
    assert m["n_entries_above_raw"] == 1
    assert m["excess_counts"] == 5
    assert m["n_cells_with_excess"] == 1
    assert m["n_genes_with_excess"] == 1


def test_new_nonzeros(raw):
    sub = raw[:20].copy()
    X = sub.X.toarray()
    X[0, 0] = 0
    raw2 = _adata(X.copy(), list(sub.obs_names), list(sub.var_names))
    X[0, 0] = 2
    out = _adata(X, list(sub.obs_names), list(sub.var_names))
    assert compute_metrics(raw2, out)["n_new_nonzeros"] == 1


def test_unknown_cells_and_float_output(raw):
    sub = raw[:5].copy()
    out = _adata(sub.X.toarray().astype(float), list(sub.obs_names[:4]) + ["other"], list(sub.var_names))
    m = compute_metrics(raw, out)
    assert m["n_cells_not_in_raw"] == 1
    assert m["is_integer"] is True


# ---------------------------------------------------------------- barnyard

@pytest.fixture
def barnyard():
    # 4 genes: h1, h2 (human), m1, m2 (mouse); 5 cells
    X = np.array([
        [50, 30, 4, 6],    # c0 human
        [40, 40, 2, 8],    # c1 human
        [5, 5, 60, 30],    # c2 mouse
        [30, 30, 30, 30],  # c3 multiplet
        [1, 1, 1, 1],      # c4 no DEMUX label
    ])
    raw = _adata(X, [f"c{i}" for i in range(5)], ["h1", "h2", "m1", "m2"],
                 species=["human", "human", "mouse", "mouse"])
    demux = pd.DataFrame({
        "pool_id": "p", "cell_id": ["c0", "c1", "c2", "c3"],
        "species": ["human", "human", "mouse", "human"],
        "sample_id": ["p_human", "p_human", "p_mouse", None],
        "multiplet": [False, False, False, True],
    })
    return raw, demux


def test_barnyard_passthrough(barnyard):
    raw, demux = barnyard
    m = compute_metrics(raw, raw.copy(), demux=demux)
    assert m["barnyard"] is True
    assert m["n_singlets"] == 3 and m["n_multiplets_skipped"] == 1 and m["n_cells_no_label"] == 1
    assert m["human_n_cells"] == 2 and m["mouse_n_cells"] == 1
    assert m["human_xspecies_frac_raw"] == pytest.approx(20 / 180)
    assert m["mouse_xspecies_frac_raw"] == pytest.approx(10 / 100)
    assert m["all_xspecies_frac_raw"] == pytest.approx(30 / 280)
    assert m["all_xspecies_frac_out"] == pytest.approx(30 / 280)
    assert m["all_xspecies_removed"] == pytest.approx(0)
    assert m["all_own_removed"] == pytest.approx(0)


def test_barnyard_perfect_removal(barnyard):
    raw, demux = barnyard
    X = raw.X.toarray().astype(float)
    X[0:2, 2:4] = 0      # human cells: mouse genes gone
    X[2, 0:2] = 0        # mouse cell: human genes gone
    X[0, 0] -= 5         # a little own-species loss
    out = _adata(X, list(raw.obs_names), list(raw.var_names))  # output var has no species column
    m = compute_metrics(raw, out, demux=demux)
    assert m["all_xspecies_removed"] == pytest.approx(1)
    assert m["all_xspecies_frac_out"] == pytest.approx(0)
    assert m["human_own_removed"] == pytest.approx(5 / 160)
    assert m["mouse_own_removed"] == pytest.approx(0)
    assert m["all_cell_xspecies_frac_out_median"] == pytest.approx(0)


def test_barnyard_needs_both_species(barnyard):
    raw, demux = barnyard
    raw.var["species"] = "mouse"
    m = compute_metrics(raw, raw.copy(), demux=demux)
    assert m["barnyard"] is False and m["all_xspecies_frac_raw"] is None


# ---------------------------------------------------------------- CLI

def test_cli_output_arg(tmp_path):
    base = ["--output_dir", str(tmp_path), "--name", "m", "--prepared_h5ad", "r.h5ad",
            "--filtered_h5ad", "f.h5ad", "--demultiplexed_tsv", "d.tsv"]
    args = parse_args(base + ["--output", "ambient_metrics_json=/x/y.json"])
    assert str(output_path(args)) == "/x/y.json"
    assert output_path(parse_args(base)) == tmp_path / "m_ambient_metrics.json"
