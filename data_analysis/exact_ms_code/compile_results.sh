#!/bin/bash

# compile all results
# est_procs=("lasso_glm" "lasso-SS_glm" "lasso-KF_glm" \
#     "SL_SL" "SL-SS_SL" "SPVIM_SL" "baseSL_SL" "SPVIM-RR_SL")
# algos=("1" "2")
est_procs=("lasso-LJ_glm" "lasso-BI-BL_glm" "SPVIM-RR_SL")
algos=("1")
output_dir="../results/data_analysis/"

for est_proc in ${est_procs[@]}; do
    IFS='_' read -ra select_reg_array <<< "$est_proc"
    select=${select_reg_array[0]}
    regress=${select_reg_array[1]}
    select_dir=${select//;/-}
    for alg in ${algos[@]}; do
        Rscript 06_compile_results.R --M 10 --est-type $regress --selection-type $select --algorithm $alg --output-dir $output_dir
    done
done
