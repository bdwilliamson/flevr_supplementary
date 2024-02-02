# Supplementary materials for the `flevr` paper

This repository contains code to reproduce the analyses in ["Flexible variable selection in the presence of missing data"](https://arxiv.org/abs/2202.12989) by Williamson and Huang (_International Journal of Biostatistics_, 2023). All analyses were implemented in the freely available R programming language; specifically, version 4.0.2. All analyses use the R package `flevr` version 0.0.2.

This README file provides an overview of the code available in the repository.

## Code directory

We have separated our code further into two sub-directories based on the two main objectives of the manuscript:

1. Numerical experiments to evaluate the operating characteristics of our proposed method under varying data-generating mechanisms (`sims`).
2. Developing a biomarker panel for pancreatic cancer early detection (`data_analysis`).

All analyses were performed on a Linux cluster using the Slurm batch scheduling system. The head node of the batch scheduler allows the shorthand "ml" in place of "module load". If you use a different batch scheduling system, the individual code files are flagged with the line where you can change batch variables. If you prefer to run the analyses locally, you may -- however, these analyses will then take a large amount of time.

-----

## Issues

If you encounter any bugs or have any specific questions about the analysis, please
[file an issue](https://github.com/bdwilliamson/flevr_supplementary/issues).
