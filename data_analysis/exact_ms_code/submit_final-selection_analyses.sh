#!/bin/bash

# Submit the final data analysis (for final set of selected variables)

# run the correct R script
echo -e \
    '#!/bin/bash\n Rscript 05_final_selected_vars.R ' \
    '--M ${1} --selection-type ${2} --output-dir ${3}' > run_main-select.sh
chmod u+x run_main-select.sh

# selection + regression procedures
selects=("lasso" "lasso-SS" "lasso-KF" \
    "SL" "SL-SS" "SPVIM" "SPVIM-RR")

for select in ${selects[@]}; do
    io_prefix="<path to preferred i/o directory>/output_${select}"
    save_dir="<path to preferred output directory>/data_analysis/output_${select}"

    mkdir -p $io_prefix
    io_file="$io_prefix/slurm-%A_%a.out"

    sbatch --time=7-0 -c5 --array=1-2 -e $io_file \
        -o $io_file ./run_main-select.sh 10 $select $save_dir
done

# clean up
rm run_main-select.sh
