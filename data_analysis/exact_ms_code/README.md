# The (nearly) exact code used in the data analysis for the `flevr` paper

This repository contains code to reproduce the analyses in ["Flexible variable selection in the presence of missing data"](https://arxiv.org/abs/2202.12989) by Williamson and Huang (2022+). All analyses were implemented in the freely available R programming language; specifically, version 4.0.2. All analyses use the R package `flevr` version 0.0.2. The code will not run unless you have access to the study data, which are not publicly available. Several lines of code have been changed to maintain confidentiality for labs within the EDRN.

The R scripts contain the code necessary to do assess the performance of the variable selection procedures and to select the final set of variables:
* `clean_data.R`: clean two analysis datasets (one for each outcome)
* `00_utils.R`: functions that are useful across the analysis
* `01_impute_missing_data.R`: impute missing data
* `02_variable_selection.R`: perform variable selection based on a specified procedure
* `03_prediction_performance.R`: assess the prediction performance of a set of selected variables
* `04_run_analysis.R`: run a single replication of the cross-validated performance estimation procedure
* `05_main.R`: run many replicates of the procedure in `04_run_analysis.R`
* `05_final_selected_vars.R`: obtain the final set of selected variables
* `06_compile_results.R`: compile many sets of results
* `07_create_plots_intrinsic.R`: create plots and tables with output

The bash scripts run the analysis and compile results:
* `submit_analysis.sh`: submit an analysis for a given procedure
* `submit_analyses.sh`: submit the prediction performance analysis for all procedures
* `submit_final-selection_analyses.sh`: submit the final variable selection analysis for each procedure

Running the following code will submit all of the simulations to a Slurm scheduler:
```{bash}
chmod u+x *.sh
./submit_analyses.sh
./submit_final-selection_analyses.sh
```
