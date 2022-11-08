#!/bin/bash

# submit all sims for the ms
# first, compute the number of jobs in the job array (if I haven't passed one in)
if [ "${11}" == "" ]; then
  num_n_p_miss=$9
  njobs=`expr $2 / $3 \* $num_n_p_miss`
  arry="1-$njobs"
else
  arry=${11}
fi

if [ ${12} -eq 1 ]; then
  io_prefix="/fh/scratch/delete90/huang_y/bwillia2/$8"
else
  io_prefix="/home/bwillia2/$8"
fi
mkdir -p $io_prefix
io_file="$io_prefix/slurm-%A_%a.out"

# Takes command-line arguments
# 1: simulation name (e.g., "binomial-linear-normal-nested")
# 2: number of total reps
# 3: number of reps per job
# 4: number of boostrap reps (for stability selection)
# 5: number of MI reps
# 6: estimator ('lasso' or 'SL')
# 7: extra layer ('' or 'SS' or 'knockoffs')
# 8: use scratch for temp results
# 9: run missing-data methods or complete-case
# 10: random number seed multiplier (for if jobs fail)
if [[ ${1} =~ "correlated" ]]; then
  echo -e \
    '#!/bin/bash\n Rscript investigate_lasso_performance_intrinsic.R --sim-name $1' \
    '--nreps-total $2 --nreps-per-job $3 --b $4 --m $5 --est-type $6' \
    '--extra-layer $7 --use-scratch $8 --run-missing $9' \
    '--seed-multiplier ${10}' > call_sim_feature_select.sh
else
  echo -e \
    '#!/bin/bash\n Rscript run_sim_intrinsic_feature_select.R --sim-name $1' \
    '--nreps-total $2 --nreps-per-job $3 --b $4 --m $5 --est-type $6' \
    '--extra-layer $7 --use-scratch $8 --run-missing $9' \
    '--seed-multiplier ${10}' > call_sim_feature_select.sh
fi
chmod u+x call_sim_feature_select.sh

# run the sim
if [ ${10} -eq 0 ]; then
  sbatch --time=14-0 -c4 --array=$arry -e $io_file -o $io_file \
    ./call_sim_feature_select.sh $1 $2 $3 $4 $5 $6 $7 ${12} ${13} ${14}
else
  sbatch --qos=restart-new --partition=restart-new --time=7-0 \
    --array=$arry -e $io_file -o $io_file \
    ./call_sim_feature_select.sh $1 $2 $3 $4 $5 $6 $7 ${12} ${13} ${14}
fi
# clean up
rm call_sim_feature_select.sh
