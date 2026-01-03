# clear workspace
rm(list = ls())

# import libraries
library(ggplot2)
library(latex2exp)
library(knitr)

# run the main simulation code
clnum <- 8 # number of cores for parallel computing
times <- 200 # number of replications

n.list <- c(250, 1000, 500)
d.list <- c(2, 2, 10)
for(i in 1:3){
  n <- n.list[i]
  d <- d.list[i]
  source("./code_functions/simulation_case1.R")
  source("./code_functions/simulation_case2.R")
  source("./code_functions/simulation_case3.R")


 mianresult_list <- list(
  "case1.true.OTR"= matrix(c(case1.q.25,0,case1.q.50,0), nrow = times, ncol = 4, byrow = TRUE),
   "case1.SCL.linear" = case1.mianresult.SCL.linear,
   "case1.SCL.gaussian" = case1.mianresult.SCL.gaussian,
   "case1.wang" = case1.mianresult.wang,
   "case1.QIQ" = case1.mianresult.QIQ,
   "case2.true.OTR"= matrix(c(case2.q.25,0,case2.q.50,0), nrow = times, ncol = 4, byrow = TRUE),
   "case2.SCL.linear" = case2.mianresult.SCL.linear,
   "case2.SCL.gaussian" = case2.mianresult.SCL.gaussian,
   "case2.wang" = case2.mianresult.wang,
   "case2.QIQ" = case2.mianresult.QIQ,
  "case3.true.OTR"= matrix(c(case3.q.25,0,case3.q.50,0), nrow = times, ncol = 4, byrow = TRUE),
   "case3.SCL.linear" = case3.mianresult.SCL.linear,
   "case3.SCL.gaussian" = case3.mianresult.SCL.gaussian,
   "case3.wang" = case3.mianresult.wang,
   "case3.QIQ" = case3.mianresult.QIQ
 )
 methods <- rep(c("True OTR","SCL-Linear", "SCL-Gaussian", "Wangs method", "QIQ-learning"), 3)
 cases <- rep(c("Case 1", "Case 2", "Case 3"), each = 5)
 n_methods <- length(methods)

 result_table_compact <- data.frame(
   Method = methods,
   Case = cases,
   stringsAsFactors = FALSE
 )

 result_table_compact[["Value for tau=0.25"]] <- character(n_methods)
 result_table_compact[["MR for tau=0.25"]] <- character(n_methods)
 result_table_compact[["Value for tau=0.50"]] <- character(n_methods)
 result_table_compact[["MR for tau=0.50"]] <- character(n_methods)

 for (j in 1:n_methods) {
   method <- methods[j]
   data_matrix <- mianresult_list[[j]]
   
   means <- colMeans(data_matrix, na.rm = TRUE)
   sds <- apply(data_matrix, 2, sd, na.rm = TRUE)
   
   result_table_compact$Method[j] <- method
   
   result_table_compact[j, "Value for tau=0.25"] <- sprintf("%.3f (%.3f)", means[1], sds[1])
   result_table_compact[j, "MR for tau=0.25"] <- sprintf("%.3f (%.3f)", means[2], sds[2])
   result_table_compact[j, "Value for tau=0.50"] <- sprintf("%.3f (%.3f)", means[3], sds[3])
   result_table_compact[j, "MR for tau=0.50"] <- sprintf("%.3f (%.3f)", means[4], sds[4])
 }
  cat(kable(result_table_compact, format = "simple", caption = paste("Table S", i, ": Mean of value functions and misclassification rates for n =", n, "and d =", d)), sep = "\n", file = paste0("./output/TableS", i, ".txt"))
}
