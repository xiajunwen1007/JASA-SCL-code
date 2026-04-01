# clear workspace
rm(list = ls())

# import libraries
library(knitr)

# run the realdata analysis code
clnum <- 16 # number of cores for parallel computing
times <- 200 # number of replications

source("./code_functions/realdata_value.R")
source("./code_functions/realdata_DC.R")

# output results of Table 2 in the main paper
mianresult_list <- list(
    "SCL-Linear" = realdata.value.SCL.linear,
    "SCL-Gaussian" = realdata.value.SCL.gaussian,
    "Wangs method" = realdata.value.wang,
    "QIQ-learning" = realdata.value.QIQ
)

methods <- names(mianresult_list)
n_methods <- length(methods)

result_table_compact <- data.frame(
  Method = methods,
  stringsAsFactors = FALSE
)

result_table_compact[["Value for tau=0.25"]] <- character(n_methods)
result_table_compact[["Value for tau=0.50"]] <- character(n_methods)
result_table_compact[["Value for tau=0.75"]] <- character(n_methods)

for (j in 1:n_methods) {
  method <- methods[j]
  data_matrix <- mianresult_list[[method]]
  
  means <- colMeans(data_matrix, na.rm = TRUE)
  sds <- apply(data_matrix, 2, sd, na.rm = TRUE)
  
  result_table_compact$Method[j] <- method
  
  result_table_compact[j, "Value for tau=0.25"] <- sprintf("%.1f (%.1f)", means[1], sds[1])
  result_table_compact[j, "Value for tau=0.50"] <- sprintf("%.1f (%.1f)", means[2], sds[2])
  result_table_compact[j, "Value for tau=0.75"] <- sprintf("%.1f (%.1f)", means[3], sds[3])
}

cat(kable(result_table_compact, format = "simple", caption = "Table 2: Mean of value function for different methods on real data"), sep = "\n", file = "./output/Table2.txt")

# output results of Table 3 in the main paper
mianresult_list <- list(
    "SCL-Linear" = realdata.dc.SCL.linear,
    "SCL-Gaussian" = realdata.dc.SCL.gaussian,
    "Wangs method" = realdata.dc.wang,
    "QIQ-learning" = realdata.dc.QIQ
)

methods <- names(mianresult_list)
n_methods <- length(methods)

result_table_compact <- data.frame(
  Method = methods,
  stringsAsFactors = FALSE
)

result_table_compact[["DC for tau=0.25"]] <- character(n_methods)
result_table_compact[["DC for tau=0.50"]] <- character(n_methods)
result_table_compact[["DC for tau=0.75"]] <- character(n_methods)

for (j in 1:n_methods) {
  method <- methods[j]
  data_matrix <- mianresult_list[[method]]
  
  means <- colMeans(data_matrix, na.rm = TRUE)
  sds <- apply(data_matrix, 2, sd, na.rm = TRUE)
  
  result_table_compact$Method[j] <- method
  
  result_table_compact[j, "DC for tau=0.25"] <- sprintf("%.3f (%.3f)", means[1], sds[1])
  result_table_compact[j, "DC for tau=0.50"] <- sprintf("%.3f (%.3f)", means[2], sds[2])
  result_table_compact[j, "DC for tau=0.75"] <- sprintf("%.3f (%.3f)", means[3], sds[3])
}

cat(kable(result_table_compact, format = "simple", caption = "Table 3: Mean of DC for different methods on real data"), sep = "\n", file = "./output/Table3.txt")