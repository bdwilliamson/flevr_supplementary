#!/bin/bash

# Submit the data analysis

# selection + regression procedures
est_procs=("lasso_glm" "lasso-SS_glm" "lasso-KF_glm" \
    "SL_SL" "SL-SS_SL" "SPVIM_SL" "SPVIM-RR_SL" "baseSL_SL")
algos=("1" "2")

# args are:
# 01: the variable selection procedure to use (e.g., "lasso")
# 02: the estimation procedure to use (e.g., "glm")
# 03: the algorithm to use (e.g., "1")
# 04: the total number of reps
# 05: the number of reps per job
# 06: the number of imputations
# 07: the number of CV folds
# 08: the output destination
# 09: the i/o prefix
# 10: the list of jobs
# 11: covariates + biomarkers?
for est_proc in ${est_procs[@]}; do
    IFS='_' read -ra select_reg_array <<< "$est_proc"
    select=${select_reg_array[0]}
    regress=${select_reg_array[1]}
    select_dir=${select//;/-}
    if [ $select == "SL-SS" ] || [ $select == "SPVIM" ] || [ $select == "SPVIM-RR" ]; then
        nreps_per_job=5
    else
        nreps_per_job=10
    fi
    for alg in ${algos[@]}; do
        io_prefix="mi_predictiveness/output_${select_dir}_${regress}_${alg}"
        # edit the next line to your preferred directory for saving output
        save_dir="<path to preferred output directory>/data_analysis/output_${select_dir}_${regress}_${alg}"
        if [ $select == "SPVIM-RR" ] && [ $alg == "2" ]; then
            # do nothing
            :
        else
            ./submit_analysis.sh $select $regress $alg 100 $nreps_per_job 10 5 $save_dir $io_prefix "" 0
        fi
    done
done
