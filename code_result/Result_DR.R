# clear workspace
rm(list = ls())

# import libraries
library(ggplot2)
library(latex2exp)
library(knitr)

# set sample size and dimension, and run the main simulation code
clnum <- 16 # number of cores for parallel computing
times <- 200 # number of replications
n <- 500
d <- 2

source("./code_functions/simulation_case1_dr.R")
source("./code_functions/simulation_case2_dr.R")
source("./code_functions/simulation_case3_dr.R")

# generate Tables S11, S12, S13 in the supplementary material, respectively.
for (i in 1:3){  
  mianresult_list <- list(
    "True OTR" = matrix(c(eval(parse(text = paste0("case", i, ".q.25"))), NA, 
                             eval(parse(text = paste0("case", i, ".q.50"))), NA), nrow = times, ncol = 4, byrow = TRUE),
    "SCL-Linear-correctcorrect" = eval(parse(text = paste0("case", i, ".SCL.linear.correctcorrect"))),
    "SCL-Gaussian-correctcorrect" = eval(parse(text = paste0("case", i, ".SCL.gaussian.correctcorrect"))),
    "SCL-Linear-incorrect2correct" = eval(parse(text = paste0("case", i, ".SCL.linear.incorrect2correct"))),
    "SCL-Gaussian-incorrect2correct" = eval(parse(text = paste0("case", i, ".SCL.gaussian.incorrect2correct"))),
    "SCL-Linear-correctincorrect" = eval(parse(text = paste0("case", i, ".SCL.linear.correctincorrect"))),
    "SCL-Gaussian-correctincorrect" = eval(parse(text = paste0("case", i, ".SCL.gaussian.correctincorrect"))),
    "SCL-Linear-incorrect2incorrect" = eval(parse(text = paste0("case", i, ".SCL.linear.incorrect2incorrect"))),
    "SCL-Gaussian-incorrect2incorrect" = eval(parse(text = paste0("case", i, ".SCL.gaussian.incorrect2incorrect")))
  )
  methods <- c("True OTRs", rep(c("SCL-Linear", "SCL-Gaussian"), 4))
  modelspecification <- c(NA,c("Correct PS and correct SM", "Correct PS and correct SM",
                          "Incorrect PS and correct SM", "Incorrect PS and correct SM",
                          "Correct PS and incorrect SM", "Correct PS and incorrect SM",
                          "Incorrect PS and incorrect SM", "Incorrect PS and incorrect SM"))
  n_methods <- length(methods)

  result_table_compact <- data.frame(
    Method = methods,
    Model_Specification = modelspecification,
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

  cat(kable(result_table_compact, format = "simple", caption = paste0("Table S", 10 + i, ": Mean of value functions and MR under different model specifications in Case ", i, ".")), sep = "\n", 
      file = paste0("./output/TableS", 10 + i, ".txt"))
}

# generate Tables S14, S15, S16 in the supplementary material, respectively.
for (i in 1:3){  
  mianresult_list <- list(
    "True OTR" = matrix(c(eval(parse(text = paste0("case", i, ".q.25"))), NA, 
                             eval(parse(text = paste0("case", i, ".q.50"))), NA), nrow = times, ncol = 4, byrow = TRUE),
    "SCL-Linear-ps1" = eval(parse(text = paste0("case", i, ".SCL.linear.incorrect1correct"))),
    "SCL-Gaussian-ps1" = eval(parse(text = paste0("case", i, ".SCL.gaussian.incorrect1correct"))),
    "SCL-Linear-ps2" = eval(parse(text = paste0("case", i, ".SCL.linear.correctcorrect"))),
    "SCL-Gaussian-ps2" = eval(parse(text = paste0("case", i, ".SCL.gaussian.correctcorrect"))),
    "SCL-Linear-ps3" = eval(parse(text = paste0("case", i, ".SCL.linear.incorrect2correct"))),
    "SCL-Gaussian-ps3" = eval(parse(text = paste0("case", i, ".SCL.gaussian.incorrect2correct")))
  )
  methods <- c("True quantile OTR", rep(c("SCL-Linear", "SCL-Gaussian"), 3))
  PSpecification <- c(NA, "(i)", "(i)", "(ii)", "(ii)", "(iii)", "(iii)")
  n_methods <- length(methods)

  result_table_compact <- data.frame(
    Method = methods,
    Model_Specification = PSpecification,
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

  cat(kable(result_table_compact, format = "simple", caption = paste0("Table S", 13 + i, ": Mean of value functions and MR under different propensity score specifications in Case ", i, ".")), sep = "\n", file = paste0("./output/TableS", 13 + i, ".txt"))
}
