# clear workspace
rm(list = ls())

# import libraries
library(ggplot2)
library(latex2exp)
library(knitr)

# set sample size and dimension, and run the main simulation code
clnum <- 16 # number of cores for parallel computing
times <- 200 # number of replications
n <- 500 # sample size
d <- 2 # dimension of covariates
C0.list <- c(3.5, 2) # censoring rate: 20% and 40%

# construct tables with a return
cat(paste0("Table S", 19, ": Mean of value functions and MR under different model specifications for survival data\n"), file = paste0("./output/TableS", 19, ".txt"))

for (i in 1:2){
  C0 <- C0.list[i]
  
  source("./code_functions/simulation_survival.R")

  mianresult_list <- list(
    "True OTR" = matrix(c(eval(parse(text = paste0("case", 1, ".q.25"))), 0, 
                             eval(parse(text = paste0("case", 1, ".q.50"))), 0), nrow = times, ncol = 4, byrow = TRUE),
    "SCL-Linear-correctcorrect" = eval(parse(text = paste0("case", 1, ".SCL.linear.correctcorrect"))),
    "SCL-Gaussian-correctcorrect" = eval(parse(text = paste0("case", 1, ".SCL.gaussian.correctcorrect"))),
    "SCL-Linear-incorrectcorrect" = eval(parse(text = paste0("case", 1, ".SCL.linear.incorrectcorrect"))),
    "SCL-Gaussian-incorrectcorrect" = eval(parse(text = paste0("case", 1, ".SCL.gaussian.incorrectcorrect"))),
    "SCL-Linear-correctincorrect" = eval(parse(text = paste0("case", 1, ".SCL.linear.correctincorrect"))),
    "SCL-Gaussian-correctincorrect" = eval(parse(text = paste0("case", 1, ".SCL.gaussian.correctincorrect"))),
    "SCL-Linear-incorrectincorrect" = eval(parse(text = paste0("case", 1, ".SCL.linear.incorrectincorrect"))),
    "SCL-Gaussian-incorrectincorrect" = eval(parse(text = paste0("case", 1, ".SCL.gaussian.incorrectincorrect")))
  )
  methods <- c("True OTRs", rep(c("SCL-Linear", "SCL-Gaussian"), 4))
  modelspecification <- c(NA,c("Correct CM and correct SM", "Correct CM and correct SM",
                          "Incorrect CM and correct SM", "Incorrect CM and correct SM",
                          "Correct CM and incorrect SM", "Correct CM and incorrect SM",
                          "Incorrect CM and incorrect SM", "Incorrect CM and incorrect SM"))
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

  if(i==1){
    cat(kable(result_table_compact, format = "simple", caption = paste0("Censoring rates: ", 20)), sep = "\n", file = paste0("./output/TableS", 19, ".txt"), append = TRUE)
  } else {
  cat(kable(result_table_compact, format = "simple", caption = paste0("Censoring rates: ", 40)), sep = "\n", file = paste0("./output/TableS", 19, ".txt"), append = TRUE)
  }
}

# plot the true OTRs
df <- expand.grid(x1 = seq(-2, 2, length.out = 500),
                  x2 = seq(-2, 2, length.out = 500))
x <- cbind(1,df,deparse.level = 0)
names(x)[2] <- "x.1"
names(x)[3] <- "x.2"

for (i in 1:1){
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
    ggsave(paste0("./output/true_25_survival.png"), width = 7.9, height = 6)

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
    ggsave(paste0("./output/true_50_survival.png"), width = 7.9, height = 6)
}
