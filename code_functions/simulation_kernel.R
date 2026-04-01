library(MASS)
library(parallel)
library(glmnet)
library(WeightSVM)
library(quantoptr)

# source necessary functions
source("./code_functions/function_kernel.R")

### The main simulation function that generates data, fits different methods, and evaluates their performance.
### The input l denotes the index of the simulation run, the inputs-n, d, base-are used to generate data, and the inputs-base.test, x.test, opt.q.25, opt.q.50-are used to evaluate different methods. 
### The outputs include all metrics for competing methods, including the value functions and misclassification rates.
simulation <- function(l, n, d, base, base.test, x.test, opt.q.25, opt.q.50) {
  # set the sample size and the dimension of the coefficients
  n <- n
  d <- d

  # generate data
  x <- cbind(1, array(rnorm(n * d), c(n, d)))
  ps <- exp(-0.5 - 0.5 * (x[, 2] + x[, 3])) / (1 + exp(-0.5 - 0.5 * (x[, 2] + x[, 3])))
  a <- sapply(ps, rbinom, n = 1, size = 1)

  mean.y <- expit(eval(parse(text = base)) +  a * (1 + x[, 2] + x[, 3]^2))
  y <- rnbinom(n, size = 2 * a + 1, prob = mean.y)

  data <- list(y = y, x = x, a = a)

  ### Estimate nuisance functions
  # estimation of the propensity scores
  f.a.x <- glmnet(x[, -1], a, family = binomial(link = "logit"), lambda=0)
  mean.a <- predict(f.a.x, x[, -1], type = "response")

  # estimation of the survival functions
  dataframe0 <- as.data.frame(list(y=y[a==0],x=x[a==0,-1]))
  dataframe1 <- as.data.frame(list(y=y[a==1],x=x[a==1,-1]))

  formula <- as.formula("y ~ x.1 + x.2")
  f.y.ax0 <- glm.nb(formula, data = dataframe0)
  f.y.ax1 <- glm.nb(formula, data = dataframe1)

  # successive classification learning
  kappa <- sd(y) / sqrt(n) / 6
  epsilon <- 1 / sqrt(n) / 2

  # Laplace kernel
  res.SCL.exp50.25 <- SCL.count.exponential(data, 0.25, kappa, epsilon, mean.a, f.y.ax1, f.y.ax0, kernel=0.5)
  class(res.SCL.exp50.25) <- "SCL"
  # Exponential15 kernel 
  res.SCL.exp75.25 <- SCL.count.exponential(data, 0.25, kappa, epsilon, mean.a, f.y.ax1, f.y.ax0, kernel=0.75)
  class(res.SCL.exp75.25) <- "SCL"
  # Gaussian kernel
  res.SCL.exp100.25 <- SCL.count.exponential(data, 0.25, kappa, epsilon, mean.a, f.y.ax1, f.y.ax0, kernel=1)
  class(res.SCL.exp100.25) <- "SCL"

  res.SCL.exp50.50 <- SCL.count.exponential(data, 0.50, kappa, epsilon, mean.a, f.y.ax1, f.y.ax0, kernel=0.5)
  class(res.SCL.exp50.50) <- "SCL"
  res.SCL.exp75.50 <- SCL.count.exponential(data, 0.50, kappa, epsilon, mean.a, f.y.ax1, f.y.ax0, kernel=0.75)
  class(res.SCL.exp75.50) <- "SCL"
  res.SCL.exp100.50 <- SCL.count.exponential(data, 0.50, kappa, epsilon, mean.a, f.y.ax1, f.y.ax0, kernel=1)
  class(res.SCL.exp100.50) <- "SCL"

  # successive classification learning with beta distribution used to smoothing
  res.SCL.exp50.25.beta <- SCL.count.exponential.beta(data, 0.25, kappa, epsilon, mean.a, f.y.ax1, f.y.ax0, 2, 0.5, kernel=0.5)
    class(res.SCL.exp50.25.beta) <- "SCL"
  res.SCL.exp75.25.beta <- SCL.count.exponential.beta(data, 0.25, kappa, epsilon, mean.a, f.y.ax1, f.y.ax0, 2, 0.5, kernel=0.75)
  class(res.SCL.exp75.25.beta) <- "SCL"
  res.SCL.exp100.25.beta <- SCL.count.exponential.beta(data, 0.25, kappa, epsilon, mean.a, f.y.ax1, f.y.ax0, 2, 0.5, kernel=1)
  class(res.SCL.exp100.25.beta) <- "SCL"
  
  res.SCL.exp50.50.beta <- SCL.count.exponential.beta(data, 0.50, kappa, epsilon, mean.a, f.y.ax1, f.y.ax0, 2, 0.5, kernel=0.5)
  class(res.SCL.exp50.50.beta) <- "SCL"
  res.SCL.exp75.50.beta <- SCL.count.exponential.beta(data, 0.50, kappa, epsilon, mean.a, f.y.ax1, f.y.ax0, 2, 0.5, kernel=0.75)
  class(res.SCL.exp75.50.beta) <- "SCL"
  res.SCL.exp100.50.beta <- SCL.count.exponential.beta(data, 0.50, kappa, epsilon, mean.a, f.y.ax1, f.y.ax0, 2, 0.5, kernel=1)
  class(res.SCL.exp100.50.beta) <- "SCL"

  # calculate the predicted optimal treatment for testing set
  treat.SCL.exp50.25 <- as.numeric(predict(res.SCL.exp50.25, x.test[, -1])) - 1
  treat.SCL.exp75.25 <- as.numeric(predict(res.SCL.exp75.25, x.test[, -1])) - 1
  treat.SCL.exp100.25 <- as.numeric(predict(res.SCL.exp100.25, x.test[, -1])) - 1

  treat.SCL.exp50.50 <- as.numeric(predict(res.SCL.exp50.50, x.test[, -1])) - 1
  treat.SCL.exp75.50 <- as.numeric(predict(res.SCL.exp75.50, x.test[, -1])) - 1
  treat.SCL.exp100.50 <- as.numeric(predict(res.SCL.exp100.50, x.test[, -1])) - 1

  treat.SCL.exp50.25.beta <- as.numeric(predict(res.SCL.exp50.25.beta, x.test[, -1])) - 1
  treat.SCL.exp75.25.beta <- as.numeric(predict(res.SCL.exp75.25.beta, x.test[, -1])) - 1
  treat.SCL.exp100.25.beta <- as.numeric(predict(res.SCL.exp100.25.beta, x.test[, -1])) - 1

  treat.SCL.exp50.50.beta <- as.numeric(predict(res.SCL.exp50.50.beta, x.test[, -1])) - 1
  treat.SCL.exp75.50.beta <- as.numeric(predict(res.SCL.exp75.50.beta, x.test[, -1])) - 1
  treat.SCL.exp100.50.beta <- as.numeric(predict(res.SCL.exp100.50.beta, x.test[, -1])) - 1

  # calculate the misspecification rate
  mr.SCL.exp50.25 <- 1 - sum(treat.SCL.exp50.25 == opt.q.25) / dim(x.test)[1]
  mr.SCL.exp75.25 <- 1 - sum(treat.SCL.exp75.25 == opt.q.25) / dim(x.test)[1]
  mr.SCL.exp100.25 <- 1 - sum(treat.SCL.exp100.25 == opt.q.25) / dim(x.test)[1]

  mr.SCL.exp50.50 <- 1 - sum(treat.SCL.exp50.50 == opt.q.50) / dim(x.test)[1]
  mr.SCL.exp75.50 <- 1 - sum(treat.SCL.exp75.50 == opt.q.50) / dim(x.test)[1]
  mr.SCL.exp100.50 <- 1 - sum(treat.SCL.exp100.50 == opt.q.50) / dim(x.test)[1]

  mr.SCL.exp50.25.beta <- 1 - sum(treat.SCL.exp50.25.beta == opt.q.25) / dim(x.test)[1]
  mr.SCL.exp75.25.beta <- 1 - sum(treat.SCL.exp75.25.beta == opt.q.25) / dim(x.test)[1]
  mr.SCL.exp100.25.beta <- 1 - sum(treat.SCL.exp100.25.beta == opt.q.25) / dim(x.test)[1]

  mr.SCL.exp50.50.beta <- 1 - sum(treat.SCL.exp50.50.beta == opt.q.50) / dim(x.test)[1]
  mr.SCL.exp75.50.beta <- 1 - sum(treat.SCL.exp75.50.beta == opt.q.50) / dim(x.test)[1]
  mr.SCL.exp100.50.beta <- 1 - sum(treat.SCL.exp100.50.beta == opt.q.50) / dim(x.test)[1]

  # calculate the smoothed value function constructed by uniform distribution
  v.SCL.exp50.25 <- smooth.quantile(rnbinom(100000, 2 * treat.SCL.exp50.25 + 1, expit(eval(parse(text = base.test)) + treat.SCL.exp50.25 * (1 + x.test[, 2] + x.test[, 3]^2))), 0.25)
  v.SCL.exp50.50 <- smooth.quantile(rnbinom(100000, 2 * treat.SCL.exp50.50 + 1, expit(eval(parse(text = base.test)) + treat.SCL.exp50.50 * (1 + x.test[, 2] + x.test[, 3]^2))), 0.50)

  v.SCL.exp75.25 <- smooth.quantile(rnbinom(100000, 2 * treat.SCL.exp75.25 + 1, expit(eval(parse(text = base.test)) + treat.SCL.exp75.25 * (1 + x.test[, 2] + x.test[, 3]^2))), 0.25)
  v.SCL.exp75.50 <- smooth.quantile(rnbinom(100000, 2 * treat.SCL.exp75.50 + 1, expit(eval(parse(text = base.test)) + treat.SCL.exp75.50 * (1 + x.test[, 2] + x.test[, 3]^2))), 0.50)

  v.SCL.exp100.25 <- smooth.quantile(rnbinom(100000, 2 * treat.SCL.exp100.25 + 1, expit(eval(parse(text = base.test)) + treat.SCL.exp100.25 * (1 + x.test[, 2] + x.test[, 3]^2))), 0.25)
  v.SCL.exp100.50 <- smooth.quantile(rnbinom(100000, 2 * treat.SCL.exp100.50 + 1, expit(eval(parse(text = base.test)) + treat.SCL.exp100.50 * (1 + x.test[, 2] + x.test[, 3]^2))), 0.50)

  v.SCL.exp50.25.beta <- smooth.quantile(rnbinom(100000, 2 * treat.SCL.exp50.25.beta + 1, expit(eval(parse(text = base.test)) + treat.SCL.exp50.25.beta * (1 + x.test[, 2] + x.test[, 3]^2))), 0.25)
  v.SCL.exp50.50.beta <- smooth.quantile(rnbinom(100000, 2 * treat.SCL.exp50.50.beta + 1, expit(eval(parse(text = base.test)) + treat.SCL.exp50.50.beta * (1 + x.test[, 2] + x.test[, 3]^2))), 0.50)

  v.SCL.exp75.25.beta <- smooth.quantile(rnbinom(100000, 2 * treat.SCL.exp75.25.beta + 1, expit(eval(parse(text = base.test)) + treat.SCL.exp75.25.beta * (1 + x.test[, 2] + x.test[, 3]^2))), 0.25)
  v.SCL.exp75.50.beta <- smooth.quantile(rnbinom(100000, 2 * treat.SCL.exp75.50.beta + 1, expit(eval(parse(text = base.test)) + treat.SCL.exp75.50.beta * (1 + x.test[, 2] + x.test[, 3]^2))), 0.50)
  
  v.SCL.exp100.25.beta <- smooth.quantile(rnbinom(100000, 2 * treat.SCL.exp100.25.beta + 1, expit(eval(parse(text = base.test)) + treat.SCL.exp100.25.beta * (1 + x.test[, 2] + x.test[, 3]^2))), 0.25)
  v.SCL.exp100.50.beta <- smooth.quantile(rnbinom(100000, 2 * treat.SCL.exp100.50.beta + 1, expit(eval(parse(text = base.test)) + treat.SCL.exp100.50.beta * (1 + x.test[, 2] + x.test[, 3]^2))), 0.50)

  # calculate the smoothed value function constructed by beta distribution
  betav.SCL.exp50.25 <- quantile(rnbinom(100000, 2 * treat.SCL.exp50.25 + 1, expit(eval(parse(text = base.test)) + treat.SCL.exp50.25 * (1 + x.test[, 2] + x.test[, 3]^2)))-rbeta(100000, 2, 0.5), 0.25)
  betav.SCL.exp50.50 <- quantile(rnbinom(100000, 2 * treat.SCL.exp50.50 + 1, expit(eval(parse(text = base.test)) + treat.SCL.exp50.50 * (1 + x.test[, 2] + x.test[, 3]^2)))-rbeta(100000, 2, 0.5), 0.50)

  betav.SCL.exp75.25 <- quantile(rnbinom(100000, 2 * treat.SCL.exp75.25 + 1, expit(eval(parse(text = base.test)) + treat.SCL.exp75.25 * (1 + x.test[, 2] + x.test[, 3]^2)))-rbeta(100000, 2, 0.5), 0.25)
  betav.SCL.exp75.50 <- quantile(rnbinom(100000, 2 * treat.SCL.exp75.50 + 1, expit(eval(parse(text = base.test)) + treat.SCL.exp75.50 * (1 + x.test[, 2] + x.test[, 3]^2)))-rbeta(100000, 2, 0.5), 0.50)

  betav.SCL.exp100.25 <- quantile(rnbinom(100000, 2 * treat.SCL.exp100.25 + 1, expit(eval(parse(text = base.test)) + treat.SCL.exp100.25 * (1 + x.test[, 2] + x.test[, 3]^2)))-rbeta(100000, 2, 0.5), 0.25)
  betav.SCL.exp100.50 <- quantile(rnbinom(100000, 2 * treat.SCL.exp100.50 + 1, expit(eval(parse(text = base.test)) + treat.SCL.exp100.50 * (1 + x.test[, 2] + x.test[, 3]^2)))-rbeta(100000, 2, 0.5), 0.50)

  betav.SCL.exp50.25.beta <- quantile(rnbinom(100000, 2 * treat.SCL.exp50.25.beta + 1, expit(eval(parse(text = base.test)) + treat.SCL.exp50.25.beta * (1 + x.test[, 2] + x.test[, 3]^2)))-rbeta(100000, 2, 0.5), 0.25)
  betav.SCL.exp50.50.beta <- quantile(rnbinom(100000, 2 * treat.SCL.exp50.50.beta + 1, expit(eval(parse(text = base.test)) + treat.SCL.exp50.50.beta * (1 + x.test[, 2] + x.test[, 3]^2)))-rbeta(100000, 2, 0.5), 0.50)

  betav.SCL.exp75.25.beta <- quantile(rnbinom(100000, 2 * treat.SCL.exp75.25.beta + 1, expit(eval(parse(text = base.test)) + treat.SCL.exp75.25.beta * (1 + x.test[, 2] + x.test[, 3]^2)))-rbeta(100000, 2, 0.5), 0.25)
  betav.SCL.exp75.50.beta <- quantile(rnbinom(100000, 2 * treat.SCL.exp75.50.beta + 1, expit(eval(parse(text = base.test)) + treat.SCL.exp75.50.beta * (1 + x.test[, 2] + x.test[, 3]^2)))-rbeta(100000, 2, 0.5), 0.50)

  betav.SCL.exp100.25.beta <- quantile(rnbinom(100000, 2 * treat.SCL.exp100.25.beta + 1, expit(eval(parse(text = base.test)) + treat.SCL.exp100.25.beta * (1 + x.test[, 2] + x.test[, 3]^2)))-rbeta(100000, 2, 0.5), 0.25)
  betav.SCL.exp100.50.beta <- quantile(rnbinom(100000, 2 * treat.SCL.exp100.50.beta + 1, expit(eval(parse(text = base.test)) + treat.SCL.exp100.50.beta * (1 + x.test[, 2] + x.test[, 3]^2)))-rbeta(100000, 2, 0.5), 0.50)

  return(list(
    SCL.exp50 = c(v.SCL.exp50.25, betav.SCL.exp50.25, mr.SCL.exp50.25, v.SCL.exp50.50, betav.SCL.exp50.50, mr.SCL.exp50.50),
    SCL.exp75 = c(v.SCL.exp75.25, betav.SCL.exp75.25, mr.SCL.exp75.25, v.SCL.exp75.50, betav.SCL.exp75.50, mr.SCL.exp75.50),
    SCL.exp100 = c(v.SCL.exp100.25, betav.SCL.exp100.25, mr.SCL.exp100.25, v.SCL.exp100.50, betav.SCL.exp100.50, mr.SCL.exp100.50),
    SCL.exp50.beta = c(v.SCL.exp50.25.beta, betav.SCL.exp50.25.beta, mr.SCL.exp50.25.beta, v.SCL.exp50.50.beta, betav.SCL.exp50.50.beta, mr.SCL.exp50.50.beta),
    SCL.exp75.beta = c(v.SCL.exp75.25.beta, betav.SCL.exp75.25.beta, mr.SCL.exp75.25.beta, v.SCL.exp75.50.beta, betav.SCL.exp75.50.beta, mr.SCL.exp75.50.beta),
    SCL.exp100.beta = c(v.SCL.exp100.25.beta, betav.SCL.exp100.25.beta, mr.SCL.exp100.25.beta, v.SCL.exp100.50.beta, betav.SCL.exp100.50.beta, mr.SCL.exp100.50.beta)
  ))
}

# Specifiy the data generating process. If you want to test other data generating processes, please change the base, n, err.test accordingly.
# Note that the first column of x is the intercept 1.
set.seed(0)
case3.base <- "-2.5-0.5*x[,2]+0.25*x[,3]^2"
case3.base.test <-"-2.5-0.5*x.test[,2]+0.25*x.test[,3]^2"

# produce testing set
m <- 100000 # sample size
x <- cbind(1, array(rnorm(m * d), c(m, d)))
ps <- exp(-0.5 - 0.5 * (x[, 2] + x[, 3])) / (1 + exp(-0.5 - 0.5 * (x[, 2] + x[, 3])))
a <- sapply(ps, rbinom, n = 1, size = 1)
mean.y <- expit(eval(parse(text = case3.base)) +  a * (1 + x[, 2] + x[, 3]^2))
y <- rnbinom(m, size = 2 * a + 1, prob = mean.y)

# calculate the true optimal smoothed quantile 
case3.q.25 <- uniroot(opt.q.case3, c(-1, 100), tau = 0.25)$root
case3.q.50 <- uniroot(opt.q.case3, c(-1, 100), tau = 0.50)$root

# calculate the true optimal smoothed quantile with beta distribution
case3.q.25.beta <- uniroot(opt.q.case3.beta, c(-1, 100), tau = 0.25, shape1=2, shape2=0.5)$root
case3.q.50.beta <- uniroot(opt.q.case3.beta, c(-1, 100), tau = 0.50, shape1=2, shape2=0.5)$root

# calculate the true optimal treatment for testing set
case3.opt.q.25 <- opt.reg.case3(case3.q.25)
case3.opt.q.50 <- opt.reg.case3(case3.q.50)

# parallel setup
cl <- makeCluster(getOption("cl.cores", clnum))
clusterSetRNGStream(cl, 333)
clusterEvalQ(cl, {
  library(MASS)
  library(glmnet)
  library(WeightSVM)
  library(quantoptr)
  library(quantreg)
  library(kernlab)
})
clusterExport(cl, varlist = list("SCL.count.exponential", "SCL.count.exponential.beta", "smooth.quantile", "expit", "compute_kernel_matrix", "predict.SCL"), envir = .GlobalEnv)

# parallel computing
system.time({
  res <- parLapply(cl, 1:times, simulation, n = n, d=d, base = case3.base, base.test = case3.base.test, x.test = x, opt.q.25 = case3.opt.q.25, opt.q.50 = case3.opt.q.50)
})
stopCluster(cl)

# extract results
# value funtion and misclassification rates
case3.SCL.exp50 <- array(0, dim = c(times, 6))
case3.SCL.exp75 <- array(0, dim = c(times, 6))
case3.SCL.exp100 <- array(0, dim = c(times, 6))

case3.SCL.exp50.beta <- array(0, dim = c(times, 6))
case3.SCL.exp75.beta <- array(0, dim = c(times, 6))
case3.SCL.exp100.beta <- array(0, dim = c(times, 6))

for (j in 1:times) {
  case3.SCL.exp50[j, ] <- res[[j]]$SCL.exp50
  case3.SCL.exp75[j, ] <- res[[j]]$SCL.exp75
  case3.SCL.exp100[j, ] <- res[[j]]$SCL.exp100

  case3.SCL.exp50.beta[j, ] <- res[[j]]$SCL.exp50.beta
  case3.SCL.exp75.beta[j, ] <- res[[j]]$SCL.exp75.beta
  case3.SCL.exp100.beta[j, ] <- res[[j]]$SCL.exp100.beta
}