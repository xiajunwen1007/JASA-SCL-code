# SCL: Successive classification learning for estimating quantile optimal treatment regimes

This directory contains all scripts to reproduce the simulation and real-data results. The entry point is `Run_all.R`, which sequentially sources the individual result scripts. All the tables and figures in the main paper and the supplementary material will appear in the `output/` folder after running the scripts.

## Summary of the files
- `required_packages.R`: script to install all the necessary R packages to reproduce our results
- `Run_all.R`: the main script to run all the result scripts
- `code_result/`: folder containing all the individual result scripts for simulation and real-data analysis  
  - `Result_main_500.R` produces main simulation results with sample size n=500 and covariates dimension d=2, i.e., Table 1 in the main paper, Tables S4–S8 in the supplementary material, and Figures S2–S7 in the supplementary material.  
  - `Result_main_diff_np.R` produces simulation results with different (n, d) combinations, including (250, 2), (1000, 2), and (500, 10), i.e., Tables S1–S3 in the supplementary material.  
  - `Result_inconsistency.R` produces simulation results on the inconsistency issue of Wang's method, i.e., Figure 1 in the main paper and Figure S1 in the supplementary material.  
  - `Result_alter_smoothing.R` produces simulation results with alternative smoothing techniques, i.e., Tables S9 and S10 in the supplementary material.  
  - `Result_DR.R` produces simulation results about the doubly robust property of our SCL method, i.e., Tables S11–S16 in the supplementary material.  
  - `Result_nonconvexity.R` produces simulation results for nonconvexity issue of Wang's method stated in Section S1.1 in the supplementary material.  
  - `Result_realdata.R` produces real-data analysis results on the AIDS study data, i.e., Tables 2 and 3 in the main paper.  
  - `Result_smooth_survival.R` produces plot of the survival functions and the smooth survival functions, i.e., Figure 2 in the main paper.  
  - `Result_kernels.R` produces simulation results with different kernel choices, i.e., Tables S17 and S18 in the supplementary material.  
  - `Result_survival.R` produces simulation results for survival data, i.e., Tables S19 and Figure S8 in the supplementary material.  
- `code_functions/`: folder containing all the functions used in the simulation and real-data analysis.
  - `function_main.R`: functions including our SCL method, the competing methods, and auxiliary functions
  - `function_inconsistency.R`: functions for the inconsistency issue of Wang's method
  - `function_nonconvexity.R`: functions for the nonconvexity issue of Wang's method
  - `function_kernel.R`: functions to test the robustness of our SCL method with different kernel choices
  - `function_survival.R`: functions to implement our SCL method for survival data
  - `simulation_case1.R`: functions to generate data for case 1 in the main paper and execute the competing methods
  - `simulation_case2.R`: functions to generate data for case 2 in the main paper and execute the competing methods
  - `simulation_case3.R`: functions to generate data for case 3 in the main paper and execute the competing methods
  - `simulation_case1_dr.R`: functions to test the doubly robust property of our SCL method in case 1
  - `simulation_case2_dr.R`: functions to test the doubly robust property of our SCL method in case 2
  - `simulation_case3_dr.R`: functions to test the doubly robust property of our SCL method in case 3
  - `simulation_inconsistency.R`: functions on the inconsistency issue of Wang's method
  - `simulation_nonconvexity.R`: functions on the nonconvexity issue of Wang's method
  - `simulation_kernel.R`: functions to generate data and test the robustness of our SCL method with different kernel choices
  - `simulation_survival.R`: functions to generate survival data and test our SCL method for survival data
  - `simulation_alter_smoothing.R`: functions to test the alternative smoothing techniques for discrete outcomes for our SCL method
  - `realdata_value.R`: functions to clean the AIDS study data (real data) and implement the cross-validated to evaluate competing methods in value function 
  - `realdata_DC.R`: functions to clean the AIDS study data (real data) and implement the cross-validated to evaluate competing methods in decision concordance (DC)
- `output/`: folder to save all the tables and figures in the main paper and the supplementary material

## Prerequisites
### R version
R version 4.3.1.
### Multi-core parallelization on a single machine
Number of cores used: 16.
### Necessary packages
The necessary packages to reproduce our results are as follows. You can run the script `required_packages.R` to intall all these packages.
- cowplot_1.1.3
- doParallel_1.0.17
- dplyr_1.1.4
- foreach_1.5.2
- ggplot2_3.5.1
- ggtext_0.1.2
- glmnet_4.1.8
- grid_4.3.1
- knitr_1.43
- latex2exp_0.9.6
- MASS_7.3.60
- mpath_0.4.2.25
- parallel_4.3.1
- quantoptr_0.1.3
- quantreg_5.97
- speff2trial_1.0.5
- tidyr_1.3.1
- WeightSVM_1.7.11
- survival_3.5.5

## How to replicate our simulation and real-data results
Execute:
	```r
	source("./Run_all.R")
	```

This will run:

  - `code_result/Result_main_500.R`
  - `code_result/Result_main_diff_np.R`
  - `code_result/Result_alter_smoothing.R`
  - `code_result/Result_inconsistency.R`
  - `code_result/Result_DR.R`
  - `code_result/Result_nonconvexity.R`
  - `code_result/Result_realdata.R`
  - `code_result/Result_smooth_survival.R`
  - `code_result/Result_kernels.R`
  - `code_result/Result_survival.R`

All the tables and figures will be saved in the `output/` folder.