#!/bin/bash

# -----------------------------------------------------------
# create plots from sim results
# -----------------------------------------------------------
# accepts the following command-line arguments
# 1: sim name (e.g., "binomial-linear-normal-nested")
# 2: number of bootstrap replicates
# 3: number of imputes
# 4: number of total jobs
# 5: should we put all missing data in? 0/1
# 6: should we restrict to lasso-based and intrinsic selection?
# 7: should we run complete-data only?
sim_name=$1
b=$2
m=$3
nreps_total=$4

output_prefix="compile_output/"
mkdir -p ${output_prefix}

Rscript create_plots.R --sim-name "$sim_name" --b "$b" --m "$m" --nreps-total "$nreps_total" --all-miss $5 --restrict-ests $6 --complete-only $7 > "${output_prefix}${1}.out" 2>&1
