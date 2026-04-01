# source and run all result scripts
source("./code_result/Result_main_500.R")
source("./code_result/Result_alter_smoothing.R")
source("./code_result/Result_inconsistency.R")
source("./code_result/Result_DR.R")
source("./code_result/Result_nonconvexity.R")
source("./code_result/Result_realdata.R")
source("./code_result/Result_smooth_survival.R")
source("./code_result/Result_kernel.R")
source("./code_result/Result_survival.R")

# Reproducing the results of Tables S1--S3 in the Supplementary Material is computationally intensive. Therefore, to reduce execution time, we have commented out the following line. To replicate these results, please uncomment it.
#source("./code_result/Result_main_diff_np.R")
