# clear workspace
rm(list = ls())

# import libraries
library(knitr)

# set sample size and dimension, and run the main simulation code
clnum <- 16 # number of cores for parallel computing
times <- 200 # number of replications
n <- 500
d <- 2
source("./code_functions/simulation_kernel.R")

# output results of Table S9 in the supplementary material
mianresult_list <- list(
    "True OTR"= matrix(c(case3.q.25, case3.q.25.beta, 0), nrow = times, ncol = 3, byrow = TRUE),
    "SCL-Laplace" = case3.SCL.exp50[,1:3],
    "SCL-Laplace.beta" = case3.SCL.exp50.beta[,1:3],
    "SCL-Exponential75" = case3.SCL.exp75[,1:3],
    "SCL-Exponential75.beta" = case3.SCL.exp75.beta[,1:3],
    "SCL-Gaussian" = case3.SCL.exp100[,1:3],
    "SCL-Gaussian.beta" = case3.SCL.exp100.beta[,1:3]
)

methods <- names(mianresult_list)
n_methods <- length(methods)

result_table_compact <- data.frame(
  Method = methods,
  stringsAsFactors = FALSE
)

result_table_compact[["Value1 for tau=0.25"]] <- character(n_methods)
result_table_compact[["Value2 for tau=0.25"]] <- character(n_methods)
result_table_compact[["MR for tau=0.25"]] <- character(n_methods)

for (j in 1:n_methods) {
  method <- methods[j]
  data_matrix <- mianresult_list[[method]]
  
  means <- colMeans(data_matrix, na.rm = TRUE)
  sds <- apply(data_matrix, 2, sd, na.rm = TRUE)
  
  result_table_compact$Method[j] <- method
  
  result_table_compact[j, "Value1 for tau=0.25"] <- sprintf("%.3f (%.3f)", means[1], sds[1])
  result_table_compact[j, "Value2 for tau=0.25"] <- sprintf("%.3f (%.3f)", means[2], sds[2])
  result_table_compact[j, "MR for tau=0.25"] <- sprintf("%.3f (%.3f)", means[3], sds[3])
}

cat(kable(result_table_compact, format = "simple", caption = "Table S17: Mean of value functions and MR for different kernels in Case 3 when quantile level tau=0.25."), sep = "\n", file = "./output/TableS17.txt")

# output results of Table S10 in the supplementary material
mianresult_list <- list(
    "True OTR"= matrix(c(case3.q.50, case3.q.50.beta, 0), nrow = times, ncol = 3, byrow = TRUE),
    "SCL-Laplace" = case3.SCL.exp50[,4:6],
    "SCL-Laplace.beta" = case3.SCL.exp50.beta[,4:6],
    "SCL-Exponential75" = case3.SCL.exp75[,4:6],
    "SCL-Exponential75.beta" = case3.SCL.exp75.beta[,4:6],
    "SCL-Gaussian" = case3.SCL.exp100[,4:6],
    "SCL-Gaussian.beta" = case3.SCL.exp100.beta[,4:6]
)

methods <- names(mianresult_list)
n_methods <- length(methods)

result_table_compact <- data.frame(
  Method = methods,
  stringsAsFactors = FALSE
)

result_table_compact[["Value1 for tau=0.50"]] <- character(n_methods)
result_table_compact[["Value2 for tau=0.50"]] <- character(n_methods)
result_table_compact[["MR for tau=0.50"]] <- character(n_methods)

for (j in 1:n_methods) {
  method <- methods[j]
  data_matrix <- mianresult_list[[method]]
  
  means <- colMeans(data_matrix, na.rm = TRUE)
  sds <- apply(data_matrix, 2, sd, na.rm = TRUE)
  
  result_table_compact$Method[j] <- method
  
  result_table_compact[j, "Value1 for tau=0.50"] <- sprintf("%.3f (%.3f)", means[1], sds[1])
  result_table_compact[j, "Value2 for tau=0.50"] <- sprintf("%.3f (%.3f)", means[2], sds[2])
  result_table_compact[j, "MR for tau=0.50"] <- sprintf("%.3f (%.3f)", means[3], sds[3])
}

cat(kable(result_table_compact, format = "simple", caption = "Table S18: Mean of value functions and MR for different kernels in Case 3 when quantile level tau=0.50."), sep = "\n", file = "./output/TableS18.txt")
