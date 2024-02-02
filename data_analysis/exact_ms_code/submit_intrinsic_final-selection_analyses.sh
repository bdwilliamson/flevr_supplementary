#!/bin/bash

# Submit the final data analysis (for final set of selected variables)
ml fhR/4.0.2-foss-2019b
ml jbigkit

# run the correct R script
echo -e \
    '#!/bin/bash\n Rscript 05_final_selected_vars.R ' \
    '--M ${1} --selection-type ${2} --output-dir ${3}' > run_main-select.sh
chmod u+x run_main-select.sh

# selection + regression procedures
selects=("lasso-LJ" "lasso-BI-BL" "SPVIM-RR")

for select in ${selects[@]}; do
    io_prefix="/fh/scratch/delete90/huang_y/bwillia2/mi_predictiveness/output_${select}"
    save_dir="/fh/fast/huang_y/bwillia2/mi_predictiveness/data_analysis/output_${select}"

    mkdir -p $io_prefix
    io_file="$io_prefix/slurm-%A_%a.out"

    sbatch --time=14-0 -c5 --array=1-4 -e $io_file \
        -o $io_file ./run_main-select.sh 10 $select $save_dir
done

# clean up
rm run_main-select.sh
