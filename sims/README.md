# Running the numerical experiments for the `flevr` paper

This repository contains code to reproduce the analyses in ["Flexible variable selection in the presence of missing data"](https://arxiv.org/abs/2202.12989) by Williamson and Huang (2022+). All analyses were implemented in the freely available R programming language; specifically, version 4.0.2. All analyses use the R package `flevr` version 0.0.2.

## Properties under various data-generating mechanisms

The following code will replicate the results in Section 4 of the main manuscript (Figures 1 and 2) and Section 2 of the Supplementary Material.

The simulation uses the following files:
* `submit_all_intrinsic_ms_sims.sh`: Submit all simulations for this section.
* `submit_sim_intrinsic_feature_select.sh`: Batch submission of a group of jobs to a Slurm scheduler.
* `run_sim_intrinsic_feature_select.R`: the main R script for this simulation, corresponding to Scenarios 1 and 3--5. Runs the simulation `nreps_per_job` times for a specified set of parameters.
* `investigate_lasso_performance_intrinsic.R`: the second main R script for this simulation, corresponding to Scenarios 2 and 6--8. Runs the simulation `nreps_per_job` times for a specified set of parameters.
* `do_one_intrinsic.R`: Runs the simulation a single time for a specified set of parameters.
* `gen_data.R`: Generate a dataset.
* `get_true_performance.R`: Get prediction performance of a selected set of variables on independent data.
* `utils.R`: Utility functions.
* `00_utils.R`: Other utility functions, shared with `data_analysis`.
* `02_variable_selection.R`: Useful functions for Long + Johnson variable selection, shared with `data_analysis`.

Running the following code will submit all of the simulations to a Slurm scheduler:
```{bash}
chmod u+x *.sh
./submit_all_intrinsic_ms_sims.sh
```
If you aren't running on a Slurm scheduler, make sure to edit the appropriate lines (flagged in each file). You can run this code locally, but it will take some time.

Once you have run the simulations and copied the results to directory `sim_output`, the following code reproduces all plots and tables:
* `create_all_plots.sh`: creates all plots for the main manuscript and supplement
* `create_plots.sh`: run `create_plots.R` for a given scenario
* `create_plots.R`: load in the results and create figures for variable selection and prediction performance

```{bash}
./create_all_plots.sh
```