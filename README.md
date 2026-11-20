# ambient-metrics

Metrics for the ambient RNA removal benchmark ([ambient-plan](https://github.com/omni-scrna/ambient-plan)).

An [omnibenchmark](https://omnibenchmark.org) module for the stage `AMBIENT_M`, which runs after `AMBIENT`, once per pool and ambient-removal method. It compares the method's output with the raw counts of the same cells and writes one flat JSON of numbers.

## Entrypoint `default` (`run.py`)

| flag | meaning |
|---|---|
| `--rawdata_h5ad` | raw droplets of one pool (h5ad, X = counts; `var['species']` for mixed-species pools) |
| `--filtered_h5ad` | output of one ambient-removal method for that pool |
| `--demultiplexed_tsv` | DEMUX output (`cell_id`, `species`, `multiplet`, ...) |
| `--output_dir`, `--name` | set by omnibenchmark (`--output ambient_metrics_json=PATH` with api 0.8) |
| `--seed` (0) | seed for the cell sample of the alignment check |

Output: `<name>_ambient_metrics.json`.

Cells and genes are matched **by name**, never by position. The `no_cleaning` method is the control: it must give 0 removed and 0 added counts.

### Count inflation (all datasets)
Ambient removal should only subtract counts.

| key | meaning |
|---|---|
| `n_cells_not_in_raw`, `n_genes_not_in_raw`, `is_integer`, `n_negative` | validity |
| `frac_counts_removed` | 1 − total out / total raw |
| `n_entries_above_raw`, `excess_counts`, `frac_excess_of_out`, `n_new_nonzeros` | counts the method **added** |
| `frac_cells_with_excess`, `cell_excess_frac_median` / `_p90` / `_max` | how the added counts spread over cells |
| `align_cell_corr_matched` / `_shifted`, `align_gene_corr_matched` / `_shifted` | name-independent alignment check (matched ≫ shifted) |

### Barnyard (mixed-species pools, e.g. HGMM)
In a human cell, counts on mouse genes can only be ambient RNA (and vice versa). Cell species come from DEMUX (run on the raw counts), gene species from the raw matrix. Multiplets and cells without a DEMUX label are left out. Keys exist for `all`, `human` and `mouse` cells:

| key | meaning | ideal |
|---|---|---|
| `<g>_xspecies_frac_raw` / `_out` | cross-species share of the counts, before / after | out → 0 |
| `<g>_xspecies_removed` | share of cross-species counts removed | 1 |
| `<g>_own_removed` | share of own-species counts removed | small (own-species ambient RNA exists, so not 0) |
| `<g>_cell_xspecies_frac_raw_median` / `_out_median` | per-cell median of the cross-species share | |

Pools without both species give `barnyard: false` and `None` for these keys.

## Environment and tests
The plan's `py_metrics` env (`envs/py_metrics.yml`: anndata ≥ 0.11, which reads the anndataR string columns of the raw files).
Tests: `conda env create -f environment.yml`, then `python -m pytest tests`.

`src/common/cli.py` is copied from [boilerplate](https://github.com/omni-scrna/boilerplate) (`scripts/pull.py`); don't edit it by hand. The stage schema `src/common/schema/AMBIENT_M.json` is written here because ambient-plan has no `schema/` folder yet.
