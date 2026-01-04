# SCL: Successive classification learning for estimating quantile optimal treatment regimes

This directory contains all scripts to reproduce the simulation and real-data results. The entry point is `Run_all.R`, which sets the working directory and sequentially sources the individual result scripts. All the tables and figures in the main paper and the supplementary material will appear in the `output/` folder after running the scripts.

## Summary of the files
- `required_packages.R`: script to install all the necessary R packages to reproduce our results
- `Run_all.R`: the main script to run all the result scripts
- `code_result/`: folder containing all the individual result scripts for simulation and real-data analysis
  - `Result_main_500.R`: main simulation results with n=500 and d=2
  - `Result_main_diff_np.R`: main simulation results with different (n, d) combinations, including (250, 2), (1000, 2), and (500, 10)
  - `Result_inconsistency.R`: simulation results on the inconsistency issue of Wang's method
  - `Result_alter_smoothing.R`: simulation results with alternative smoothing techniques
  - `Result_DR.R`: simulation results about the doubly robust property of our SCL method
  - `Result_nonconvexity.R`: simulation results for nonconvexity issue of Wang's method
  - `Result_realdata.R`: real-data analysis results on the AIDS study data
  - `Result_smooth_survival.R`: plot of the survival functions and the smooth survival functions.
- `code_functions/`: folder containing all the functions used in the simulation and real-data.
  - `function_main.R`: functions including our SCL method, the competing methods, and auxiliary functions
  - `simulation_case1.R`: functions to generate data for case 1 in the main paper and execute the competing methods
  - `simulation_case2.R`: functions to generate data for case 2 in the main paper and execute the competing methods
  - `simulation_case3.R`: functions to generate data for case 3 in the main paper and execute the competing methods
  - `simulation_case1_dr.R`: functions to test the doubly robust property of our SCL method in case 1
  - `simulation_case2_dr.R`: functions to test the doubly robust property of our SCL method in case 2
  - `simulation_case3_dr.R`: functions to test the doubly robust property of our SCL method in case 3
  - `simulation_inconsistency.R`: functions on the inconsistency issue of Wang's method
  - `simulation_nonconvexity.R`: functions on the nonconvexity issue of Wang's method
  - `realdata_value.R`: functions to clean the AIDS study data (real data) and implement the cross-validated to evaluate competing methods in value function 
  - `realdata_RI.R`: functions to clean the AIDS study data (real data) and implement the cross-validated to evaluate competing methods in rand index (RI)
- `code_functions/`: folder containing all the functions used in the result scripts
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

## How to replicate our simulation and real-data results
1. Change the directory `Documents/quantile_optimal_regime/submission/` in the script `Run_all.R` to your local path.
2. Execute:
	```r
	source("./Run_all_simulation.R")
	```
	This will set the working directory and run:
	- `code_result/Result_main_500.R`
	- `code_result/Result_main_diff_np.R`
	- `code_result/Result_alter_smoothing.R`
	- `code_result/Result_inconsistency.R`
	- `code_result/Result_DR.R`
	- `code_result/Result_nonconvexity.R`
	- `code_result/Result_realdata.R`
	- `code_result/Result_smooth_survival.R`
3. All the tables and figures will be saved in the `output/` folder.