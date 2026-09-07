# Manuscript result map

The table below maps each reported simulation result to its generating scripts,
source data, and manuscript-ready output. Table and figure numbers follow the
current main manuscript and Supplementary Material.

| Manuscript result | Design | Simulation command | Reported source data | Manuscript-ready output |
|---|---|---|---|---|
| Main Table 1 | Gaussian errors, aligned anomalies, IID predictors | `Rscript run_detection.R main` | `results/reported/detection_main/` | `outputs/reported/detection_main_iid.tex` |
| Main Table 2 | Gaussian errors, aligned anomalies, AR predictors | `Rscript run_detection.R main` | `results/reported/detection_main/` | `outputs/reported/detection_main_ar.tex` |
| Main Table 3 | Gaussian-error common-slope estimation, IID and AR predictors | `Rscript run_estimation.R` | `results/reported/estimation/` | `outputs/reported/estimation_normal.tex` |
| Supplementary Table S1 | Mixture-normal errors, aligned anomalies, IID predictors | `Rscript run_detection.R mixnorm` | `results/reported/detection_mixnorm/` | `outputs/reported/detection_mixnorm_iid.tex` |
| Supplementary Table S2 | Mixture-normal errors, aligned anomalies, AR predictors | `Rscript run_detection.R mixnorm` | `results/reported/detection_mixnorm/` | `outputs/reported/detection_mixnorm_ar.tex` |
| Supplementary Table S3 | AR(1) errors, aligned anomalies, IID predictors | `Rscript run_detection.R ar1` | `results/reported/detection_ar1/` | `outputs/reported/detection_ar1_iid.tex` |
| Supplementary Table S4 | AR(1) errors, aligned anomalies, AR predictors | `Rscript run_detection.R ar1` | `results/reported/detection_ar1/` | `outputs/reported/detection_ar1_ar.tex` |
| Supplementary Table S5 | Gaussian errors, random anomaly directions, IID predictors | `Rscript run_detection.R random` | `results/reported/detection_random/` | `outputs/reported/detection_random_iid.tex` |
| Supplementary Table S6 | Gaussian errors, random anomaly directions, AR predictors | `Rscript run_detection.R random` | `results/reported/detection_random/` | `outputs/reported/detection_random_ar.tex` |
| Supplementary Table S7 | Mixture-normal common-slope estimation | `Rscript run_estimation.R` | `results/reported/estimation/` | `outputs/reported/estimation_mixnorm.tex` |
| Supplementary Table S8 | AR(1)-error common-slope estimation | `Rscript run_estimation.R` | `results/reported/estimation/` | `outputs/reported/estimation_ar1.tex` |
| Supplementary Table S9 | Values of `c_C` and `k_c` in the sensitivity design | `Rscript run_sensitivity.R` | `results/reported/sensitivity/tuning_grid.csv` | Typeset directly from the tuning grid in the supplement |
| Supplementary Figure S1 | Power sensitivity under `ell=0.5`, `s_n=15`, `zeta=0.2` | `Rscript run_sensitivity.R` then `Rscript make_outputs.R` | `results/reported/sensitivity/` | `outputs/reported/sensitivity_alternative_power.pdf` |
| Supplementary Figure S2 | FDR sensitivity under the same alternative | `Rscript run_sensitivity.R` then `Rscript make_outputs.R` | `results/reported/sensitivity/` | `outputs/reported/sensitivity_alternative_fdr.pdf` |
| Supplementary Figure S3 | FDR sensitivity under the global null | `Rscript run_sensitivity.R` then `Rscript make_outputs.R` | `results/reported/sensitivity/` | `outputs/reported/sensitivity_global_null_fdr.pdf` |

## Detection-family arguments

`run_detection.R` accepts one argument:

- `main`: Gaussian errors and aligned anomaly directions.
- `mixnorm`: mixture-normal errors and aligned anomaly directions.
- `ar1`: stationary AR(1) errors with coefficient 0.5 and aligned directions.
- `random`: Gaussian errors and independent unit-sphere anomaly directions.

For `random`, the global-null rows are copied from `main` replication by
replication. This is intentional: under the global null all departure vectors
are zero, so anomaly direction is undefined and the two data-generating
processes coincide.

## Output construction

`make_outputs.R` reads the summary CSV files, requires 500 successful
replications in every reported cell, and creates all LaTeX tables and
sensitivity figures. `audit_results.R` independently reads the replication-level
and summary CSV files and checks the manuscript design before output creation.

