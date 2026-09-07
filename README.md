# Simulation code for FDR-controlled slope-anomaly detection

This repository reproduces the simulation results in the manuscript
*FDR-Controlled Detection of Anomalous Regression Slopes in Repeatedly Observed Systems*.
It contains the data-generating mechanisms, the proposed procedure, seven
comparison procedures, Monte Carlo drivers, manuscript-output builders, and the
replication-level results used in the paper.

The repository covers the simulation study only. It does not contain the Divvy
application or any real-data files.

## Quick verification

The reported 500-replication results are included under `results/reported/`.
From the repository root, audit those files and rebuild every simulation table
and sensitivity figure with:

```bash
Rscript run_all.R report
```

The regenerated files are written to `outputs/reproduced/`. This command checks
the design dimensions, tuning values, replication counts, random-direction
global-null reuse, numerical summaries, and exact LaTeX table contents.

To run a small end-to-end computational test without changing the reported
results:

```bash
Rscript run_all.R smoke
```

Smoke-test files are written to `results/smoke/` and are ignored by Git.

## Full Monte Carlo reproduction

To regenerate all experiments with 500 replications:

```bash
Rscript run_all.R full
```

Full results are written to `results/recomputed/`, and the corresponding tables
and figures are written to `outputs/reproduced/`. The scripts save each design
cell separately and resume from completed cells. A full run is computationally
intensive; runtime depends strongly on the available cores.

The number of workers can be controlled with:

```bash
DET_CORES=4 EST_CORES=4 SENS_CORES=4 Rscript run_all.R full
```

On Windows PowerShell, set the variables before calling R:

```powershell
$env:DET_CORES = "4"
$env:EST_CORES = "4"
$env:SENS_CORES = "4"
Rscript run_all.R full
```

## Simulation design

- `n = 200`, `T = 100`, and `K = m = 4`.
- All reported experiments use 500 Monte Carlo replications.
- The proposed procedure uses `C0 = c_C T^(-1/3)` with the default `c_C = 1`
  and `k_c = 0.60`.
- The anomaly counts are `s_n = ceiling(n^ell)` for
  `ell in {0.3, 0.5, 0.7}`, giving 5, 15, and 41 anomalous units.
- Detection uses `zeta in {0, 0.1, 0.2, 0.3}`; `zeta = 0` is a separately
  generated global null with every departure vector equal to zero.
- Estimation uses `ell = 0.7` and `zeta in {0.1, 0.3, 0.5}`.
- The reported predictor processes are IID and AR. No MA predictor design is
  included in this repository.
- The primary errors are Gaussian. The supplement additionally considers
  mixture-normal and AR(1) errors with `phi = 0.5`.
- The primary anomalies have aligned directions. The supplement additionally
  considers independent random directions on the unit sphere.
- Unit intercepts are generated independently from `U[-1,1]` without a
  zero-sum normalization.

Further design details are documented in `MANUSCRIPT_RESULTS.md` and directly
in the simulation drivers.

## Repository files

| File | Purpose |
|---|---|
| `run_all.R` | Unified `report`, `smoke`, and `full` entry point. |
| `simulation_core.R` | Data generation, clean-set screening, unit-level testing, BH selection, benchmark methods, and performance evaluation. |
| `study_utils.R` | Shared family specifications, parallel-worker handling, and summary utilities. |
| `run_detection.R` | Runs Gaussian, mixture-normal, AR(1)-error, or random-direction detection experiments. |
| `run_estimation.R` | Runs common-slope RMSE and MAE experiments after method-specific anomaly removal. |
| `run_sensitivity.R` | Runs the factorial sensitivity analysis for `c_C` and `k_c`. |
| `make_outputs.R` | Builds the LaTeX tables, sensitivity figures, and compact result digests. |
| `audit_results.R` | Checks all reported result dimensions, design constants, and numerical consistency conditions. |
| `validate_implementation.R` | Deterministic checks of anomaly counts, intercept generation, random directions, AR(1) errors, screening, Power, and FDR. |
| `verify_reported_outputs.R` | Compares regenerated manuscript tables and digests with the checked-in reported versions. |
| `MANUSCRIPT_RESULTS.md` | Maps every script and output to the corresponding main-text or supplementary table/figure. |
| `RESULTS_CODEBOOK.md` | Defines the columns in the replication-level and summary result files. |
| `SHA256SUMS` | SHA-256 checksums for the publishable repository files. |
| `results/reported/` | Replication-level CSV files and summaries underlying the published tables and figures. |
| `outputs/reported/` | Exact LaTeX tables and figure files used by the current manuscript. |

## Running individual experiments

The four detection families are:

```bash
Rscript run_detection.R main
Rscript run_detection.R mixnorm
Rscript run_detection.R ar1
Rscript run_detection.R random
```

The remaining experiments are:

```bash
Rscript run_estimation.R
Rscript run_sensitivity.R
```

By default, these commands write to `results/recomputed/`. Each driver accepts
environment variables for smaller diagnostic runs; see the header of the
corresponding R file.

## Software requirements

The simulation code uses only R base and recommended packages, including
`parallel`, `stats`, `graphics`, and `grDevices`. No package installation step
is required in a standard R distribution. The reported results were generated
with R 4.6.0 on Windows 11; the full session information is retained in each
reported-results directory.

## Reproducibility notes

- Every Monte Carlo replication has an explicit recorded seed.
- Random-direction experiments reuse the Gaussian global-null results exactly,
  because anomaly direction is undefined under the global null.
- The IID global-null sensitivity experiment reuses the panels underlying main
  Table 1 and verifies the default tuning result replication by replication.
- The simulated panels are regenerated from their seeds and are not stored.
  Replication-level performance outputs are retained as CSV files.
- Parallel execution does not change the seed assigned to any replication.


