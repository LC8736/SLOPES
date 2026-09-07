# Results codebook

The files under `results/reported/` contain the 500-replication results used in
the manuscript. Missing `Power` under the global null is intentional because
there are no anomalous units in that setting.

## Detection replication-level files

Each `detection_raw.csv` contains one row for each method, design cell, and
replication.

| Column | Definition |
|---|---|
| `method` | Detection procedure. |
| `selected`, `retained` | Numbers of selected and retained units. |
| `fallback` | Indicator that all units were selected and the estimation fallback was used. |
| `TP`, `FP` | Numbers of true and false discoveries. |
| `Power` | `TP / s_n`; undefined under the global null. |
| `FDR` | Realized false discovery proportion, `FP / max(selected, 1)`. Its Monte Carlo mean is the reported FDR. |
| `RMSE`, `MAE` | Common-slope estimation errors after method-specific removal; retained here for diagnostics. |
| `reference_purity` | Clean fraction of a conformal reference set, when applicable. |
| `reference_clean_recall` | Fraction of clean reference-panel units included, when applicable. |
| `scenario`, `family` | Scenario and detection-family identifiers. |
| `zeta` | Signal upper-bound multiplier; zero denotes the global null. |
| `ell`, `s_n` | Sparsity exponent and realized number of anomalous units. |
| `c0`, `k_c` | Clean-set screening parameters used by the proposed procedure. |
| `replication`, `seed` | Monte Carlo replication number and deterministic random seed. |
| `predictor` | `iid` or `ar`. |
| `error` | `normal`, `mixnorm`, or `ar1`. |
| `direction` | `aligned` or `spherical`. |
| `error_phi` | AR(1) error coefficient when applicable. |

Each `detection_summary.csv` reports Monte Carlo means and standard deviations
of Power, FDR, and the number selected, together with replication and failure
counts.

## Estimation files

`estimation_raw.csv` contains one row for each method, error design, predictor
design, signal level, and replication. `RMSE` and `MAE` are computed from the
common-slope estimator after removing the units selected by that method. The
`Raw` row uses all units without anomaly removal.

`estimation_summary.csv` contains the corresponding Monte Carlo means and
standard deviations. All estimation experiments use `ell = 0.7` and `s_n = 41`.

## Sensitivity files

`sensitivity_raw.csv` contains one row for each generated panel, predictor,
signal setting, and `(c_C, k_c)` pair. `C0` is the realized value
`c_C T^(-1/3)`. The two settings are:

- `global-null`: `zeta = 0` and `s_n = 0`;
- `alternative`: `ell = 0.5`, `s_n = 15`, and `zeta = 0.2`.

`sensitivity_summary.csv` gives Monte Carlo means and standard deviations.
`tuning_grid.csv` records the 20 parameter pairs used in Supplementary Table S9.

## Auxiliary files

- `design.csv` records the simulation cells and their checkpoint identifiers.
- `vsom_cutoffs.csv` records the predictor-specific VSOM bootstrap cutoffs.
- `sessionInfo.txt` records the R environment used for the reported run.

