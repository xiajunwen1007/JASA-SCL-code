# clear workspace
rm(list = ls())

# import libraries
library(ggplot2)
library(latex2exp)
library(dplyr)
library(tidyr)
library(cowplot)
library(ggtext)
library(grid)
library(knitr)

# set sample size and dimension, and run the main simulation code
clnum <- 16 # number of cores for parallel computing
times <- 200 # number of replications
source("./code_functions/simulation_inconsistency.R")

# prepare data for the indicator of whether the global maximizer is achieved
inconsistency_list <- list(
  "tau0.40" = hatQ.wang.40,
  "tau0.50" = hatQ.wang.50,
  "tau0.60" = hatQ.wang.60
)

methods <- c("tau=0.40", "tau=0.50", "tau=0.60")
n_methods <- length(methods)

result_table_compact <- data.frame(
  Quantile_level = methods,
  stringsAsFactors = FALSE
)

result_table_compact[["n=30"]] <- character(n_methods)
result_table_compact[["n=100"]] <- character(n_methods)
result_table_compact[["n=500"]] <- character(n_methods)
result_table_compact[["n=1000"]] <- character(n_methods)

for (j in 1:n_methods) {
  method <- methods[j]
  data_matrix <- inconsistency_list[[j]]
  
  means <- colMeans(data_matrix, na.rm = TRUE)
  sds <- apply(data_matrix, 2, sd, na.rm = TRUE)
  
  result_table_compact$Method[j] <- method
  
  result_table_compact[j, "n=30"] <- sprintf("%.3f", means[1])
  result_table_compact[j, "n=100"] <- sprintf("%.3f", means[2])
  result_table_compact[j, "n=500"] <- sprintf("%.3f", means[3])
  result_table_compact[j, "n=1000"] <- sprintf("%.3f", means[4])
}

cat(kable(result_table_compact, format = "simple", caption = "Estimated value \\( \\widehat{Q}(Y(\\hat{d})) \\) for Wang's methods under different sample sizes and tau values."), sep = "\n", file = "./output/Result_inconsistency.txt")


# prepare data for plotting tau=0.40
df <- data.frame(
  SampleSize = rep(c(30, 100, 500, 1000), times = 2),
  Method = rep(c("Wang", "SCL-Linear"), each = 4),
  Error = 1-c(colMeans(v.wang.40),colMeans(v.SCL.linear.40))
)

# create plot for MSE when tau=0.40
p1 <- ggplot(df, aes(x = SampleSize, y = Error, color = Method, group = Method)) +
  geom_line(aes(linetype = Method), size = 1) +
  geom_point(size = 3) +
  labs(x = "Sample Size", y = "MSE of quantile values") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_blank(),
    legend.position = "bottom",
    text = element_text(size = 18)
  ) +
  scale_x_continuous(breaks = c(0,100,500, 1000),labels = c("0", "100", "500", "1000")) +
  scale_color_manual(values = c("Wang" = "blue", "Smoothed Wang" = "red", "SCL-Linear" = "black")) +
  scale_linetype_manual(values = c("Wang" = "solid", "Smoothed Wang"="dashed", "SCL-Linear" = "dotted")) +
  guides(color = guide_legend(override.aes = list(shape = NA), keywidth = 5))+
  ggtitle(TeX("$\\tau=0.40$"))


# prepare data for plotting tau=0.50
df <- data.frame(
  SampleSize = rep(c(30, 100, 500, 1000), times = 2),
  Method = rep(c("Wang", "SCL-Linear"), each = 4),
  Error = 1-c(colMeans(v.wang.50),colMeans(v.SCL.linear.50))
)

# create plot for MSE when tau=0.50
p2 <- ggplot(df, aes(x = SampleSize, y = Error, color = Method, group = Method)) +
  geom_line(aes(linetype = Method), size = 1) +
  geom_point(size = 3) +
  labs(x = "Sample Size", y = "MSE of quantile values") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_blank(),
    legend.position = "bottom",
    text = element_text(size = 18)
  ) +
  scale_x_continuous(breaks = c(0, 100, 500, 1000)) +
  scale_color_manual(values = c("Wang" = "blue", "Smoothed Wang" = "red", "SCL-Linear" = "black")) +
  scale_linetype_manual(values = c("Wang" = "solid", "Smoothed Wang"="dashed", "SCL-Linear" = "dotted")) +
  guides(color = guide_legend(override.aes = list(shape = NA), keywidth = 5))+
  ggtitle(TeX("$\\tau=0.50$"))

# prepare data for plotting tau=0.60
df <- data.frame(
  SampleSize = rep(c(30, 100, 500, 1000), times = 2),
  Method = rep(c("Wang", "SCL-Linear"), each = 4),
  Error = 1-c(colMeans(v.wang.60),colMeans(v.SCL.linear.60))
)

# create plot for MSE when tau=0.60
p3 <- ggplot(df, aes(x = SampleSize, y = Error, color = Method, group = Method)) +
  geom_line(aes(linetype = Method), size = 1) +
  geom_point(size = 3) +
  labs(x = "Sample Size", y = "MSE of quantile values") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_blank(),
    legend.position = "bottom",
    text = element_text(size = 18)
  ) +
  scale_x_continuous(breaks = c(0, 100, 500, 1000)) +
  scale_color_manual(values = c("Wang" = "blue", "Smoothed Wang" = "red", "SCL-Linear" = "black")) +
  scale_linetype_manual(values = c("Wang" = "solid", "Smoothed Wang"="dashed", "SCL-Linear" = "dotted")) +
  guides(color = guide_legend(override.aes = list(shape = NA), keywidth = 5))+
  ggtitle(TeX("$\\tau=0.60$"))

# combine plots with shared legend
get_legend <- function(myggplot) {
  tmp <- ggplotGrob(myggplot)
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend <- get_legend(p1 + theme(legend.position = "bottom"))

p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")
p3 <- p3 + theme(legend.position = "none")

combined_plot <- plot_grid(p1, p2, p3, ncol = 3, align = "v")
final_plot <- plot_grid(combined_plot, legend, ncol = 1, rel_heights = c(1, 0.1))

# save the final plot 
ggsave("./output/inconsistency.png", width = 14, height = 7, dpi = 300)