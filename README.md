# Simulation code for FDR-controlled slope-anomaly detection

This repository reproduces the simulation results in the manuscript
*FDR-Controlled Detection of Anomalous Regression Slopes in Repeatedly Observed Systems*.
It contains the data-generating mechanisms, the proposed procedure, seven
comparison procedures, Monte Carlo drivers, manuscript-output builders, and the
replication-level results used in the paper.





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




## Software requirements

The simulation code uses only R base and recommended packages, including
`parallel`, `stats`, `graphics`, and `grDevices`. No package installation step
is required in a standard R distribution. The reported results were generated
with R 4.6.0 on Windows 11; the full session information is retained in each
reported-results directory.


## Real data
dizzy_real_data_code.zip contains real data and code to reproduce the results in Real Data Analysis.


