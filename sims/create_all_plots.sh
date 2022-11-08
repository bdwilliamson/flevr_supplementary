#!/bin/bash

# create all plots from sims

# -----------------------------------------------------------
# create plots for the main manuscript
# -----------------------------------------------------------
main_sim_names=("binomial-probit-linear-normal-nested" "nonlinear-normal-correlated")
for sim_name in ${main_sim_names[@]}; do
  ./create_plots.sh $sim_name 100 10 1000 0 0 0
done
# -----------------------------------------------------------
# create plots for the supplement
# -----------------------------------------------------------
supp_sim_names=("linear-normal-correlated" "linear-normal-uncorrelated" \
                "nonlinear-normal-weak-uncorrelated" \
                "binomial-probit-linear-nonnormal-nested" \
                "binomial-probit-nonlinear-normal-nested" \
                "binomial-probit-nonlinear-nonnormal-nested")
# missing-data sims
for sim_name in ${supp_sim_names[@]}; do
  ./create_plots.sh $sim_name 100 10 1000 1 0 0
done
# complete-case sims
all_sims=("${main_sim_names[@]}" "${supp_sim_names[@]}")
for sim_name in ${all_sims[@]}; do
  ./create_plots.sh $sim_name 100 1 1000 1 0 1
done
