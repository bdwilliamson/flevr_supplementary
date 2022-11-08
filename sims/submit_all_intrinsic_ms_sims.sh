#!/bin/bash

# submit all the simulations for the main manuscript (and some for the supp)

# load modules
ml fhR/4.0.2-foss-2019b
ml jbigkit

# submit all simulations that involve complete-case data
est_layers_nomissing=("lasso_none" "lasso_KF" "lasso_SS" "SPVIM_none")
# settings: A, D, B, C, E, F
sim_prefix="intrinsic-binomial-probit-"
sim_suffix="-nested"
sim_names=("linear-normal" "nonlinear-normal-correlated" \
           "linear-normal-correlated" "linear-normal-uncorrelated" "nonlinear-normal-weak-uncorrelated" \
           "linear-nonnormal" "nonlinear-normal" "nonlinear-nonnormal")
use_scratch=1
run_miss=0
seed_multiplier=1
M=1
for est_layer in "${est_layers_nomissing[@]}"; do
  IFS='_' read -ra est_layer_array <<< "$est_layer"
  est=${est_layer_array[0]}
  layer=${est_layer_array[1]}
  # if not SS or SPVIM, use restart and 5 reps per job
  if [ $layer == "SS" ] || [[ $est =~ "SPVIM" ]]; then
    nreps_per_job=2
    use_restart=0
  else
    nreps_per_job=5
    use_restart=1
  fi
  # submit
  output_est=$(echo "$est_layer" | tr '[:upper:]' '[:lower:]')
  for sim_name in "${sim_names[@]}"; do
    echo "Sim: $sim_name; est: $est_layer"
    io_prefix="mi_predictiveness/output_${output_est}_${sim_name}"
    if [[ ${sim_name} =~ "correlated" ]]; then
      full_sim_name="${sim_name}"
      num_unique_settings=4
    else
      full_sim_name="${sim_prefix}${sim_name}${sim_suffix}"
      num_unique_settings=8
      io_prefix="${io_prefix}${sim_suffix}"
    fi
    ./submit_sim_intrinsic_feature_select.sh $full_sim_name 1000 $nreps_per_job 100 $M $est $layer $io_prefix $num_unique_settings $use_restart "" $use_scratch $run_miss $seed_multiplier
  done
done

# submit all simulations that involve missing data
est_layers_missing=("lasso-LJ_none" "lasso-BI-BL_none" "SPVIM-RR_none")
use_scratch=1
run_miss=1
seed_multiplier=1
M=10
for est_layer in "${est_layers_missing[@]}"; do
  IFS='_' read -ra est_layer_array <<< "$est_layer"
  est=${est_layer_array[0]}
  layer=${est_layer_array[1]}
  # if not SS or SPVIM, use restart and 5 reps per job
  if [ $layer == "SS" ] || [[ $est =~ "SPVIM" ]]; then
    nreps_per_job=2
    use_restart=0
  else
    nreps_per_job=5
    use_restart=1
  fi
  # submit
  output_est=$(echo "$est_layer" | tr '[:upper:]' '[:lower:]')
  for sim_name in "${sim_names[@]}"; do
    echo "Sim: ${sim_name}; est: $est_layer"
    io_prefix="mi_predictiveness/output_${output_est}_${sim_name}"
    if [[ ${sim_name} =~ "correlated" ]]; then
      full_sim_name="${sim_name}"
      num_unique_settings=8
    else
      full_sim_name="${sim_prefix}${sim_name}${sim_suffix}"
      num_unique_settings=16
      io_prefix="${io_prefix}${sim_suffix}"
    fi
    ./submit_sim_intrinsic_feature_select.sh $full_sim_name 1000 $nreps_per_job 100 $M $est $layer $io_prefix $num_unique_settings $use_restart "" $use_scratch $run_miss $seed_multiplier
  done
done