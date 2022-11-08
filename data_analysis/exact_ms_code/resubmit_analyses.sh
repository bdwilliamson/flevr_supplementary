#!/bin/bash

# check and resubmit failed jobs

# load modules
ml fhR/4.0.2-foss-2019b
ml jbigkit

# selection + regression procedures
# est_procs=("lasso_glm" "lasso-SS_glm" "lasso-KF_glm" \
#     "SL_SL" "SL-SS_SL" "SPVIM_SL" "baseSL_SL")
est_procs=("lasso_glm" "lasso-SS_glm" "lasso-KF_glm" \
           "SL_SL" "baseSL_SL" "SL-SS_SL" "SPVIM_SL" "SPVIM-RR_SL")
algos=("1" "2")

outcomes=("mucinous" "high_malignancy")
covariates=("0" "1")
biomarkers=("0" "1")
n_outcomes=2
n_covar_biomarker_combos=3

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
    ub=`expr 100 / $nreps_per_job`
    for alg in ${algos[@]}; do
        echo "Est: ${est_proc}; alg: ${alg}"
        io_prefix="mi_predictiveness/output_${select_dir}_${regress}_${alg}"
        save_dir="/fh/fast/huang_y/bwillia2/mi_predictiveness/data_analysis/output_${select_dir}_${regress}_${alg}"
        jobs_to_submit=""
        if [ $select == "SPVIM-RR" ] && [ $alg == "2" ]; then
            # do nothing
            :
        else
            for outcome in ${outcomes[@]}; do
                if [ $outcome == "mucinous" ]; then
                    j=1
                else
                    j=2
                fi
                for covar in ${covariates[@]}; do
                    for biomarker in ${biomarkers[@]}; do
                        if [ $covar -eq 0 ] && [ $biomarker -eq 0 ]; then
                            # do nothing
                            :
                        else
                            file_name="${save_dir}/impute_10/perf_${outcome}_${covar}_${biomarker}_id"
                            while read -r i; do
                                [[ -f "${file_name}_${i}.rds" ]] || {
                                    ((job=$i+$ub\*($j-1)+$ub\*($n_outcomes)\*(1-$covar)\*($biomarker) + $ub\*($n_outcomes)\*($n_covar_biomarker_combos - 1)\*($covar)\*($biomarker)))
                                    if [ "$jobs_to_submit" == "" ]; then
                                        jobs_to_submit="$job"
                                    else
                                        jobs_to_submit="$jobs_to_submit,$job"
                                    fi
                                }
                            done < <(seq "$ub")
                        fi
                    done
                done
            done
            # echo $jobs_to_submit
            if [ "$jobs_to_submit" == "" ]; then
                # do nothing
                echo "No jobs to submit"
            else
                ./submit_analysis.sh $select $regress $alg 100 $nreps_per_job 10 5 $save_dir $io_prefix $jobs_to_submit 0
            fi
        fi
    done
done
