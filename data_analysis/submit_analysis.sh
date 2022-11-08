#!/bin/bash

# submit a single selection + estimation + algorithm combination


# Takes 9 command-line arguments
# 01: the variable selection procedure to use (e.g., "lasso")
# 02: the estimation procedure to use (e.g., "glm")
# 03: the algorithm to use (e.g., "1")
# 04: the total number of reps
# 05: the number of reps per job
# 06: the number of imputations
# 07: the number of CV folds
# 08: the output destination
# 09: the i/o file
# 10: the list of jobs to submit
# 11: whether or not to run the biomarker + covariate analyses
# edit the next line to your preferred i/o directory
io_prefix="/fh/scratch/delete90/huang_y/bwillia2/$9"
mkdir -p $io_prefix
io_file="$io_prefix/slurm-%A_%a.out"

# run the correct R script
echo -e \
    '#!/bin/bash\n Rscript 05_main.R ' \
    '--nreps-total $4 --nreps-per-job $5 ' \
    '--K $7 --M $6 --est-type $2 ' \
    '--selection-type $1 --algorithm $3 ' \
    '--output-dir $8' > run_main.sh
chmod u+x run_main.sh

if [ "${10}" == "" ]; then
    njobs=`expr $4 / $5 \* 4`
    arry="1-$njobs"
    if [ ${11} -eq 1 ]; then
        start=`expr $njobs + 1`
        end=`expr $4 / $5 \* 2 + $njobs`
        arry="$start-$end"
    fi
else
    arry=${10}
fi

sbatch --time=7-0 -c10 --array=$arry -e $io_file \
    -o $io_file ./run_main.sh $1 $2 $3 $4 $5 $6 $7 $8
rm run_main.sh
