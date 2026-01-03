# clear workspace
rm(list = ls())

# import libraries
library(knitr)

# set sample size and dimension, and run the main simulation code
clnum <- 16 # number of cores for parallel computing
times <- 200 # number of replications

source("./code_functions/simulation_nonconvexity.R")

# output proportion of times wang's method finds the global optimum
cat("\nProportion of times Wang's method finds the global optimum:\n", mean(flag), file = "./output/Result_nonconvexity.txt")

# output value and MR of gridsearch method and wang's method
mianresult_list <- list(
    "True OTR"= matrix(c(q.50, 0), nrow = times, ncol = 2, byrow = TRUE),
    "Gridsearch method" = gridsearch.50,
    "Wang" = wang.50
)

methods <- names(mianresult_list)
n_methods <- length(methods)

result_table_compact <- data.frame(
  Method = methods,
  stringsAsFactors = FALSE
)

result_table_compact[["Value for tau=0.50"]] <- character(n_methods)
result_table_compact[["MR for tau=0.50"]] <- character(n_methods)
for (j in 1:n_methods) {
  method <- methods[j]
  data_matrix <- mianresult_list[[method]]
  
  means <- colMeans(data_matrix, na.rm = TRUE)
  sds <- apply(data_matrix, 2, sd, na.rm = TRUE)
  
  result_table_compact$Method[j] <- method
  
  result_table_compact[j, "Value for tau=0.50"] <- sprintf("%.3f (%.3f)", means[1], sds[1])
  result_table_compact[j, "MR for tau=0.50"] <- sprintf("%.3f (%.3f)", means[2], sds[2])
}

cat("\n", file = "./output/Result_nonconvexity.txt", append = TRUE)
cat(kable(result_table_compact, format = "simple", caption = "Mean of value functions and misclassification rates for tau=0.50"), sep = "\n", file = "./output/Result_nonconvexity.txt", append = TRUE)
