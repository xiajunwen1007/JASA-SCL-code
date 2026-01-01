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
source("./code_functions/simulation_case1.R")
source("./code_functions/simulation_case2.R")
source("./code_functions/simulation_case3.R")

# output results of Table 1 in the main paper
mianresult_list <- list(
  "case1.true.quantile.OTR"= matrix(c(case1.q.25, NA, case1.q.50, NA), nrow = times, ncol = 4, byrow = TRUE),
  "case1.SCL.linear" = case1.mianresult.SCL.linear,
  "case1.SCL.gaussian" = case1.mianresult.SCL.gaussian,
  "case1.wang" = case1.mianresult.wang,
  "case1.QIQ" = case1.mianresult.QIQ,
  "case2.true.quantile.OTRs"= matrix(c(case2.q.25, NA, case2.q.50, NA), nrow = times, ncol = 4, byrow = TRUE),
  "case2.SCL.linear" = case2.mianresult.SCL.linear,
  "case2.SCL.gaussian" = case2.mianresult.SCL.gaussian,
  "case2.wang" = case2.mianresult.wang,
  "case2.QIQ" = case2.mianresult.QIQ,
  "case3.true.quantile.OTRs"= matrix(c(case3.q.25, NA, case3.q.50, NA), nrow = times, ncol = 4, byrow = TRUE),
  "case3.SCL.linear" = case3.mianresult.SCL.linear,
  "case3.SCL.gaussian" = case3.mianresult.SCL.gaussian,
  "case3.wang" = case3.mianresult.wang,
  "case3.QIQ" = case3.mianresult.QIQ
)
methods <- rep(c("True OTR", "SCL-Linear", "SCL-Gaussian", "Wangs method", "QIQ-learning"), 3)
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

cat(kable(result_table_compact, format = "simple", caption = "Table 1: Mean of value functions and misclassification rates for n = 500 in three cases."), sep = "\n", file = "./output/Table1.txt")

# output results of Table S4 in the main paper
time_list <- list(
  "case1.SCL.linear" = case1.time.SCL.linear,
  "case1.SCL.gaussian" = case1.time.SCL.gaussian,
  "case1.wang" = case1.time.wang,
  "case1.QIQ" = case1.time.QIQ,
  "case2.SCL.linear" = case2.time.SCL.linear,
  "case2.SCL.gaussian" = case2.time.SCL.gaussian,
  "case2.wang" = case2.time.wang,
  "case2.QIQ" = case2.time.QIQ,
  "case3.SCL.linear" = case3.time.SCL.linear,
  "case3.SCL.gaussian" = case3.time.SCL.gaussian,
  "case3.wang" = case3.time.wang,
  "case3.QIQ" = case3.time.QIQ
)

methods <- rep(c("SCL-Linear", "SCL-Gaussian", "Wangs method", "QIQ-learning"), 3)
cases <- rep(c("Case 1", "Case 2", "Case 3"), each = 4)
n_methods <- length(methods)

result_table_compact <- data.frame(
  Method = methods,
  Case = cases,
  stringsAsFactors = FALSE
)

result_table_compact[["Time for tau=0.25"]] <- character(n_methods)
result_table_compact[["Time for tau=0.50"]] <- character(n_methods)

for (j in 1:n_methods) {
  method <- methods[j]
  data_matrix <- time_list[[j]]
  
  means <- colMeans(data_matrix, na.rm = TRUE)
  sds <- apply(data_matrix, 2, sd, na.rm = TRUE)
  
  result_table_compact$Method[j] <- method
  
  result_table_compact[j, "Time for tau=0.25"] <- sprintf("%.3f (%.3f)", means[1], sds[1])
  result_table_compact[j, "Time for tau=0.50"] <- sprintf("%.3f (%.3f)", means[2], sds[2])
}

cat(kable(result_table_compact, format = "simple", caption = "Table S4: Computation time (in seconds) for different methods in three cases."), sep = "\n", file = "./output/TableS4.txt")

# output results of Tables S5
bias_list <- list(
  "case1.SCL.linear" = case1.q.SCL.linear,
  "case1.SCL.gaussian" = case1.q.SCL.gaussian,
  "case1.wang" = case1.q.wang,
  "case1.QIQ" = case1.q.QIQ,
  "case2.SCL.linear" = case2.q.SCL.linear,
  "case2.SCL.gaussian" = case2.q.SCL.gaussian,
  "case2.wang" = case2.q.wang,
  "case2.QIQ" = case2.q.QIQ,
  "case3.SCL.linear" = case3.q.SCL.linear,
  "case3.SCL.gaussian" = case3.q.SCL.gaussian,
  "case3.wang" = case3.q.wang,
  "case3.QIQ" = case3.q.QIQ
)

methods <- rep(c("SCL-Linear", "SCL-Gaussian", "Wangs method", "QIQ-learning"), 3)
cases <- rep(c("Case 1", "Case 2", "Case 3"), each = 4)
n_methods <- length(methods)

result_table_compact <- data.frame(
  Method = methods,
  Case = cases,
  stringsAsFactors = FALSE
)

result_table_compact[["Bias for tau=0.25"]] <- character(n_methods)
result_table_compact[["Bias for tau=0.50"]] <- character(n_methods)

for (j in 1:n_methods) {
  method <- methods[j]
  data_matrix <- bias_list[[j]]
  
  means <- colMeans(data_matrix, na.rm = TRUE)
  sds <- apply(data_matrix, 2, sd, na.rm = TRUE)
  
  result_table_compact$Method[j] <- method
  
  result_table_compact[j, "Bias for tau=0.25"] <- sprintf("%.3f (%.3f)", means[1], sds[1])
  result_table_compact[j, "Bias for tau=0.50"] <- sprintf("%.3f (%.3f)", means[2], sds[2])
}

cat(kable(result_table_compact, format = "simple", caption = "Table S5: Bias of the estimated optimal value for different methods in three cases."), sep = "\n", file = "./output/TableS5.txt")

# output results of Tables S6 and S7 in the main paper
case1.vhatv.SCL.linear <- cbind(case1.mianresult.SCL.linear[,1], case1.hatv.SCL.linear[,1], case1.mianresult.SCL.linear[,3], case1.hatv.SCL.linear[,2])
case1.vhatv.SCL.gaussian <- cbind(case1.mianresult.SCL.gaussian[,1], case1.hatv.SCL.gaussian[,1], case1.mianresult.SCL.gaussian[,3], case1.hatv.SCL.gaussian[,2])
case1.vhatv.QIQ <- cbind(case1.mianresult.QIQ[,1], case1.hatv.QIQ[,1], case1.mianresult.QIQ[,3], case1.hatv.QIQ[,2])
case1.vhatv.wang <- cbind(case1.mianresult.wang[,1], case1.hatv.wang[,1], case1.mianresult.wang[,3], case1.hatv.wang[,2])
case1.vhatv.true <- cbind(case1.q.25, case1.hatv.true[,1], case1.q.50, case1.hatv.true[,2])

case2.vhatv.SCL.linear <- cbind(case2.mianresult.SCL.linear[,1], case2.hatv.SCL.linear[,1], case2.mianresult.SCL.linear[,3], case2.hatv.SCL.linear[,2])
case2.vhatv.SCL.gaussian <- cbind(case2.mianresult.SCL.gaussian[,1], case2.hatv.SCL.gaussian[,1], case2.mianresult.SCL.gaussian[,3], case2.hatv.SCL.gaussian[,2])
case2.vhatv.QIQ <- cbind(case2.mianresult.QIQ[,1], case2.hatv.QIQ[,1], case2.mianresult.QIQ[,3], case2.hatv.QIQ[,2])
case2.vhatv.wang <- cbind(case2.mianresult.wang[,1], case2.hatv.wang[,1], case2.mianresult.wang[,3], case2.hatv.wang[,2])
case2.vhatv.true <- cbind(case2.q.25, case2.hatv.true[,1], case2.q.50, case2.hatv.true[,2])

case3.vhatv.SCL.linear <- cbind(case3.mianresult.SCL.linear[,1], case3.hatv.SCL.linear[,1], case3.mianresult.SCL.linear[,3], case3.hatv.SCL.linear[,2])
case3.vhatv.SCL.gaussian <- cbind(case3.mianresult.SCL.gaussian[,1], case3.hatv.SCL.gaussian[,1], case3.mianresult.SCL.gaussian[,3], case3.hatv.SCL.gaussian[,2])
case3.vhatv.QIQ <- cbind(case3.mianresult.QIQ[,1], case3.hatv.QIQ[,1], case3.mianresult.QIQ[,3], case3.hatv.QIQ[,2])
case3.vhatv.wang <- cbind(case3.mianresult.wang[,1], case3.hatv.wang[,1], case3.mianresult.wang[,3], case3.hatv.wang[,2])
case3.vhatv.true <- cbind(case3.q.25, case3.hatv.true[,1], case3.q.50, case3.hatv.true[,2])

mianresult_list <- list(
  "case1.True" = case1.vhatv.true,
  "case1.SCL.linear" = case1.vhatv.SCL.linear,
  "case1.SCL.gaussian" = case1.vhatv.SCL.gaussian,
  "case1.wang" = case1.vhatv.wang,
  "case1.QIQ-learning" =case1.vhatv.QIQ,
  "case2.True"= case2.vhatv.true,
  "case2.SCL.linear" = case2.vhatv.SCL.linear,
  "case2.SCL.gaussian" = case2.vhatv.SCL.gaussian,
  "case2.wang" = case2.vhatv.wang,
  "case2.QIQ-learning" =case2.vhatv.QIQ,
  "case3.True"= case3.vhatv.true,
  "case3.SCL.linear" = case3.vhatv.SCL.linear,
  "case3.SCL.gaussian" = case3.vhatv.SCL.gaussian,
  "case3.wang" = case3.vhatv.wang,
  "case3.QIQ-learning" =case3.vhatv.QIQ
)

methods <- rep(c("True OTR","SCL-Linear", "SCL-Gaussian", "Wangs method", "QIQ-learning"), 3)
cases <- rep(c("Case 1", "Case 2", "Case 3"), each = 5)
n_methods <- length(methods)

result_table_compact <- data.frame(
  Method = methods,
  Case = cases,
  stringsAsFactors = FALSE
)

result_table_compact[["True value for tau=0.25"]] <- character(n_methods)
result_table_compact[["Estimated value for tau=0.25"]] <- character(n_methods)

for (j in 1:n_methods) {
  method <- methods[j]
  data_matrix <- mianresult_list[[j]]
  
  means <- colMeans(data_matrix, na.rm = TRUE)
  sds <- apply(data_matrix, 2, sd, na.rm = TRUE)
  
  result_table_compact$Method[j] <- method
  
  result_table_compact[j, "True value for tau=0.25"] <- sprintf("%.3f (%.3f)", means[1], sds[1])
  result_table_compact[j, "Estimated value for tau=0.25"] <- sprintf("%.3f (%.3f)", means[2], sds[2])
}

cat(kable(result_table_compact, format = "simple", caption = "Table S6: True and estimated values with tau=0.25 for different methods in three cases."), sep = "\n", file = "./output/TableS6.txt")

# output results of Table S7 in the main paper
result_table_compact <- data.frame(
  Method = methods,
  Case = cases,
  stringsAsFactors = FALSE
)

result_table_compact[["True value for tau=0.50"]] <- character(n_methods)
result_table_compact[["Estimated value for tau=0.50"]] <- character(n_methods)

for (j in 1:n_methods) {
  method <- methods[j]
  data_matrix <- mianresult_list[[j]]
  
  means <- colMeans(data_matrix, na.rm = TRUE)
  sds <- apply(data_matrix, 2, sd, na.rm = TRUE)
  
  result_table_compact$Method[j] <- method
  
  result_table_compact[j, "True value for tau=0.50"] <- sprintf("%.3f (%.3f)", means[3], sds[3])
  result_table_compact[j, "Estimated value for tau=0.50"] <- sprintf("%.3f (%.3f)", means[4], sds[4])
}

cat(kable(result_table_compact, format = "simple", caption = "Table S7: True and estimated values with tau=0.50 for different methods in three cases."), sep = "\n", file = "./output/TableS7.txt")

# output table s8
mianresult_list <- list(
  "case1.true.quantile.OTR"= matrix(c(case1.q.25, NA, case1.q.50, NA), nrow = times, ncol = 4, byrow = TRUE),
  "case1.SCL.linear.soft" = case1.mianresult.SCL.linear,
  "case1.SCL.linear.hard" = case1.mianresult.SCL.linear.hard,
  "case1.SCL.gaussian.soft" = case1.mianresult.SCL.gaussian,
  "case1.SCL.gaussian.hard" = case1.mianresult.SCL.gaussian.hard,
  "case2.true.quantile.OTR"= matrix(c(case2.q.25, NA, case2.q.50, NA), nrow = times, ncol = 4, byrow = TRUE),
  "case2.SCL.linear.soft" = case2.mianresult.SCL.linear,
  "case2.SCL.linear.hard" = case2.mianresult.SCL.linear.hard,
  "case2.SCL.gaussian.soft" = case2.mianresult.SCL.gaussian,
  "case2.SCL.gaussian.hard" = case2.mianresult.SCL.gaussian.hard,
  "case3.true.quantile.OTR"= matrix(c(case3.q.25, NA, case3.q.50, NA), nrow = times, ncol = 4, byrow = TRUE),
  "case3.SCL.linear.soft" = case3.mianresult.SCL.linear,
  "case3.SCL.linear.hard" = case3.mianresult.SCL.linear.hard,
  "case3.SCL.gaussian.soft" = case3.mianresult.SCL.gaussian,
  "case3.SCL.gaussian.hard" = case3.mianresult.SCL.gaussian.hard
)
methods <- rep(c("True OTR", "SCL-Linear-soft", "SCL-Linear-hard", "SCL-Gaussian-soft", "SCL-Gaussian-hard"), 3)
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

cat(kable(result_table_compact, format = "simple", caption = "Table S8: Simulation results for SCL methods with soft and hard regimes."), sep = "\n", file = "./output/TableS8.txt")

# plot the results for Figure S2-S8 
df <- expand.grid(x1 = seq(-2, 2, length.out = 500),
                  x2 = seq(-2, 2, length.out = 500))
x <- cbind(1,df,deparse.level = 0)
names(x)[2] <- "x.1"
names(x)[3] <- "x.2"

for (i in 1:3){
  df$xgrid.treat.SCL.linear.25 <- colMeans(eval(parse(text = paste0("case", i, ".xgrid.treat.SCL.linear.25"))))
  df$xgrid.treat.SCL.gaussian.25 <- colMeans(eval(parse(text = paste0("case", i, ".xgrid.treat.SCL.gaussian.25"))))
  df$xgrid.treat.SCL.linear.50 <- colMeans(eval(parse(text = paste0("case", i, ".xgrid.treat.SCL.linear.50"))))
  df$xgrid.treat.SCL.gaussian.50 <- colMeans(eval(parse(text = paste0("case", i, ".xgrid.treat.SCL.gaussian.50"))))
  df$xgrid.treat.QIQ.25 <- colMeans(eval(parse(text = paste0("case", i, ".xgrid.treat.QIQ.25"))))
  df$xgrid.treat.QIQ.50 <- colMeans(eval(parse(text = paste0("case", i, ".xgrid.treat.QIQ.50"))))
  df$xgrid.treat.wang.25 <- colMeans(eval(parse(text = paste0("case", i, ".xgrid.treat.wang.25"))))
  df$xgrid.treat.wang.50 <- colMeans(eval(parse(text = paste0("case", i, ".xgrid.treat.wang.50"))))
  df$opt.q.25 <- eval(parse(text = paste0("opt.reg.case", i, "(case", i, ".q.25)")))
  df$opt.q.50 <- eval(parse(text = paste0("opt.reg.case", i, "(case", i, ".q.50)")))

  ggplot(df, aes(x1, x2)) +
    geom_raster(mapping = aes(fill = factor(opt.q.25))) +
    xlab(TeX("$X_1$")) +
    ylab(TeX("$X_2$")) +
    theme_minimal() +
    theme(
      legend.position = "right", 
      text = element_text(size = 22),
      legend.title = element_text(vjust = 3)
    ) +
    scale_fill_manual(
      values = c("1" = "#2A60D3", "0" = "#73C2CD"),
      breaks = c("1", "0"),
      name = "Quantile\nOptimal\nTreatment",
      labels = c("1" = "1", "0" = "0")
    )
    ggsave(paste0("./output/Case", i, "_true_25.png"), width = 7.9, height = 6)

    ggplot(df, aes(x1, x2)) +
    geom_raster(mapping = aes(fill = factor(opt.q.50))) +
    xlab(TeX("$X_1$")) +
    ylab(TeX("$X_2$")) +
    theme_minimal() +
    theme(
      legend.position = "right", 
      text = element_text(size = 22),
      legend.title = element_text(vjust = 3)
    ) +
    scale_fill_manual(
      values = c("1" = "#2A60D3", "0" = "#73C2CD"),
      name = "Quantile\nOptimal\nTreatment",
      breaks = c("1", "0"),
      labels = c("1" = "1", "0" = "0")
    )
    ggsave(paste0("./output/Case", i, "_true_50.png"), width = 7.9, height = 6)

    ggplot(df, aes(x1, x2)) +
    geom_raster(mapping=aes(fill= xgrid.treat.SCL.linear.25)) +
    xlab(TeX("$X_1$")) +
    ylab(TeX("$X_2$")) +
    theme_minimal()+
    theme(legend.position = "right", text = element_text(size = 22), legend.title = element_text(vjust = 3))+
  scale_fill_gradientn(
      colours = colorRampPalette(c("#73C2CD","#2A60D3"))(20),
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      name = "AQOT"
    )+
    geom_contour(mapping=aes(z=opt.q.25), breaks = 0.5, color = "white")

    ggsave(paste0("./output/Case", i, "_SCL.linear_25.png"), width = 7.2, height = 6)

    ggplot(df, aes(x1, x2)) +
    geom_raster(mapping=aes(fill= xgrid.treat.SCL.gaussian.25)) +
    xlab(TeX("$X_1$")) +
    ylab(TeX("$X_2$")) +
    theme_minimal()+  
    theme(legend.position = "right", text = element_text(size = 22), legend.title = element_text(vjust = 3))+
  scale_fill_gradientn(
      colours = colorRampPalette(c("#73C2CD","#2A60D3"))(20),
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      name = "AQOT"
    )+
    geom_contour(mapping=aes(z=opt.q.25), breaks = 0.5, color = "white")

    ggsave(paste0("./output/Case", i, "_SCL.gaussian_25.png"), width = 7.2, height = 6)

    ggplot(df, aes(x1, x2)) +
    geom_raster(mapping=aes(fill= xgrid.treat.SCL.linear.50)) +
    xlab(TeX("$X_1$")) +
    ylab(TeX("$X_2$")) +
    theme_minimal()+
    theme(legend.position = "right", text = element_text(size = 22), legend.title = element_text(vjust = 3))+
  scale_fill_gradientn(
      colours = colorRampPalette(c("#73C2CD","#2A60D3"))(20),
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      name = "AQOT"
    )+
    geom_contour(mapping=aes(z=opt.q.50), breaks = 0.5, color = "white")

    ggsave(paste0("./output/Case", i, "_SCL.linear_50.png"), width = 7.2, height = 6)

    ggplot(df, aes(x1, x2)) +
    geom_raster(mapping=aes(fill= xgrid.treat.SCL.gaussian.50)) +
    xlab(TeX("$X_1$")) +
    ylab(TeX("$X_2$")) +  
    theme_minimal()+
    theme(legend.position = "right", text = element_text(size = 22), legend.title = element_text(vjust = 3))+
  scale_fill_gradientn(
      colours = colorRampPalette(c("#73C2CD","#2A60D3"))(20),
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      name = "AQOT"
    )+
    geom_contour(mapping=aes(z=opt.q.50), breaks = 0.5, color = "white")

    ggsave(paste0("./output/Case", i, "_SCL.gaussian_50.png"), width = 7.2, height = 6)

    ggplot(df, aes(x1, x2)) +
    geom_raster(mapping=aes(fill= xgrid.treat.QIQ.25)) +
    xlab(TeX("$X_1$")) +
    ylab(TeX("$X_2$")) +
    theme_minimal()+
    theme(legend.position = "right", text = element_text(size = 22), legend.title = element_text(vjust = 3))+
  scale_fill_gradientn(
      colours = colorRampPalette(c("#73C2CD","#2A60D3"))(20),
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      name = "AQOT"
    )+
    geom_contour(mapping=aes(z=opt.q.25), breaks = 0.5, color = "white")

    ggsave(paste0("./output/Case", i, "_QIQlearning_25.png"), width = 7.2, height = 6)

    ggplot(df, aes(x1, x2)) +
    geom_raster(mapping=aes(fill= xgrid.treat.QIQ.50)) +
    xlab(TeX("$X_1$")) +
    ylab(TeX("$X_2$")) +
    theme_minimal()+  
    theme(legend.position = "right", text = element_text(size = 22), legend.title = element_text(vjust = 3))+
  scale_fill_gradientn(
      colours = colorRampPalette(c("#73C2CD","#2A60D3"))(20),
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      name = "AQOT"
    )+
    geom_contour(mapping=aes(z=opt.q.50), breaks = 0.5, color = "white")

    ggsave(paste0("./output/Case", i, "_QIQlearning_50.png"), width = 7.2, height = 6)

    ggplot(df, aes(x1, x2)) +
    geom_raster(mapping=aes(fill= xgrid.treat.wang.25)) +
    xlab(TeX("$X_1$")) +
    ylab(TeX("$X_2$")) +
    theme_minimal()+
    theme(legend.position = "right", text = element_text(size = 22), legend.title = element_text(vjust = 3))+
  scale_fill_gradientn(
      colours = colorRampPalette(c("#73C2CD","#2A60D3"))(20),
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      name = "AQOT"
    )+
    geom_contour(mapping=aes(z=opt.q.25), breaks = 0.5, color = "white")

    ggsave(paste0("./output/Case", i, "_Wangs.method_25.png"), width = 7.2, height = 6)

    ggplot(df, aes(x1, x2)) +
    geom_raster(mapping=aes(fill= xgrid.treat.wang.50)) +
    xlab(TeX("$X_1$")) +
    ylab(TeX("$X_2$")) +
    theme_minimal()+
    theme(legend.position = "right", text = element_text(size = 22), legend.title = element_text(vjust = 3))+
  scale_fill_gradientn(
      colours = colorRampPalette(c("#73C2CD","#2A60D3"))(20),
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      name = "AQOT"
    )+
    geom_contour(mapping=aes(z=opt.q.50), breaks = 0.5, color = "white")

    ggsave(paste0("./output/Case", i, "_Wangs.method_50.png"), width = 7.2, height = 6)
}
