library(MASS)
library(parallel)
library(glmnet)
library(WeightSVM)
library(quantoptr)

# source necessary functions
source("./code_functions/function_main.R")

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

  res.SCL.linear.25 <- SCL.count(data, 0.25, kappa, epsilon, mean.a, f.y.ax1, f.y.ax0, kernel="linear")
  res.SCL.gaussian.25 <- SCL.count(data, 0.25, kappa, epsilon, mean.a, f.y.ax1, f.y.ax0, kernel="gaussian")
  res.SCL.linear.50 <- SCL.count(data, 0.50, kappa, epsilon, mean.a, f.y.ax1, f.y.ax0, kernel="linear")
  res.SCL.gaussian.50 <- SCL.count(data, 0.50, kappa, epsilon, mean.a, f.y.ax1, f.y.ax0, kernel="gaussian")

  # successive classification learning with beta distribution used to smoothing
  res.SCL.linear.25.beta <- SCL.count.beta(data, 0.25, kappa, epsilon, mean.a, f.y.ax1, f.y.ax0, 2, 0.5, kernel="linear")
  res.SCL.gaussian.25.beta <- SCL.count.beta(data, 0.25, kappa, epsilon, mean.a, f.y.ax1, f.y.ax0, 2, 0.5, kernel="gaussian")
  res.SCL.linear.50.beta <- SCL.count.beta(data, 0.50, kappa, epsilon, mean.a, f.y.ax1, f.y.ax0, 2, 0.5, kernel="linear")
  res.SCL.gaussian.50.beta <- SCL.count.beta(data, 0.50, kappa, epsilon, mean.a, f.y.ax1, f.y.ax0, 2, 0.5, kernel="gaussian")

  # extract results
  treat.SCL.linear.25 <- as.numeric(predict(res.SCL.linear.25$opt.reg, x.test[, -1])) - 1
  treat.SCL.gaussian.25 <- as.numeric(predict(res.SCL.gaussian.25$opt.reg, x.test[, -1])) - 1
  treat.SCL.linear.50 <- as.numeric(predict(res.SCL.linear.50$opt.reg, x.test[, -1])) - 1
  treat.SCL.gaussian.50 <- as.numeric(predict(res.SCL.gaussian.50$opt.reg, x.test[, -1])) - 1

  treat.SCL.linear.25.beta <- as.numeric(predict(res.SCL.linear.25.beta$opt.reg, x.test[, -1])) - 1
  treat.SCL.gaussian.25.beta <- as.numeric(predict(res.SCL.gaussian.25.beta$opt.reg, x.test[, -1])) - 1
  treat.SCL.linear.50.beta <- as.numeric(predict(res.SCL.linear.50.beta$opt.reg, x.test[, -1])) - 1
  treat.SCL.gaussian.50.beta <- as.numeric(predict(res.SCL.gaussian.50.beta$opt.reg, x.test[, -1])) - 1

  # calculate the misspecification rate
  mr.SCL.linear.25 <- 1 - sum(treat.SCL.linear.25 == opt.q.25) / dim(x.test)[1]
  mr.SCL.gaussian.25 <- 1 - sum(treat.SCL.gaussian.25 == opt.q.25) / dim(x.test)[1]
  mr.SCL.linear.50 <- 1 - sum(treat.SCL.linear.50 == opt.q.50) / dim(x.test)[1]
  mr.SCL.gaussian.50 <- 1 - sum(treat.SCL.gaussian.50 == opt.q.50) / dim(x.test)[1]
  
  mr.SCL.linear.25.beta <- 1 - sum(treat.SCL.linear.25.beta == opt.q.25) / dim(x.test)[1]
  mr.SCL.gaussian.25.beta <- 1 - sum(treat.SCL.gaussian.25.beta == opt.q.25) / dim(x.test)[1]
  mr.SCL.linear.50.beta <- 1 - sum(treat.SCL.linear.50.beta == opt.q.50) / dim(x.test)[1]
  mr.SCL.gaussian.50.beta <- 1 - sum(treat.SCL.gaussian.50.beta == opt.q.50) / dim(x.test)[1]
  
  # calculate the smoothed value function constructed by uniform distribution
  v.SCL.linear.25 <- smooth.quantile(rnbinom(100000, 2 * treat.SCL.linear.25 + 1, expit(eval(parse(text = base.test)) + treat.SCL.linear.25 * (1 + x.test[, 2] + x.test[, 3]^2))), 0.25)
  v.SCL.linear.50 <- smooth.quantile(rnbinom(100000, 2 * treat.SCL.linear.50 + 1, expit(eval(parse(text = base.test)) + treat.SCL.linear.50 * (1 + x.test[, 2] + x.test[, 3]^2))), 0.50)

  v.SCL.gaussian.25 <- smooth.quantile(rnbinom(100000, 2 * treat.SCL.gaussian.25 + 1, expit(eval(parse(text = base.test)) + treat.SCL.gaussian.25 * (1 + x.test[, 2] + x.test[, 3]^2))), 0.25)
  v.SCL.gaussian.50 <- smooth.quantile(rnbinom(100000, 2 * treat.SCL.gaussian.50 + 1, expit(eval(parse(text = base.test)) + treat.SCL.gaussian.50 * (1 + x.test[, 2] + x.test[, 3]^2))), 0.50)

  v.SCL.linear.25.beta <- smooth.quantile(rnbinom(100000, 2 * treat.SCL.linear.25.beta + 1, expit(eval(parse(text = base.test)) + treat.SCL.linear.25.beta * (1 + x.test[, 2] + x.test[, 3]^2))), 0.25)
  v.SCL.linear.50.beta <- smooth.quantile(rnbinom(100000, 2 * treat.SCL.linear.50.beta + 1, expit(eval(parse(text = base.test)) + treat.SCL.linear.50.beta * (1 + x.test[, 2] + x.test[, 3]^2))), 0.50)

  v.SCL.gaussian.25.beta <- smooth.quantile(rnbinom(100000, 2 * treat.SCL.gaussian.25.beta + 1, expit(eval(parse(text = base.test)) + treat.SCL.gaussian.25.beta * (1 + x.test[, 2] + x.test[, 3]^2))), 0.25)
  v.SCL.gaussian.50.beta <- smooth.quantile(rnbinom(100000, 2 * treat.SCL.gaussian.50.beta + 1, expit(eval(parse(text = base.test)) + treat.SCL.gaussian.50.beta * (1 + x.test[, 2] + x.test[, 3]^2))), 0.50)

  # calculate the smoothed value function constructed by beta distribution
  betav.SCL.linear.25 <- quantile(rnbinom(100000, 2 * treat.SCL.linear.25 + 1, expit(eval(parse(text = base.test)) + treat.SCL.linear.25 * (1 + x.test[, 2] + x.test[, 3]^2)))-rbeta(100000, 2, 0.5), 0.25)
  betav.SCL.linear.50 <- quantile(rnbinom(100000, 2 * treat.SCL.linear.50 + 1, expit(eval(parse(text = base.test)) + treat.SCL.linear.50 * (1 + x.test[, 2] + x.test[, 3]^2)))-rbeta(100000, 2, 0.5), 0.50)

  betav.SCL.gaussian.25 <- quantile(rnbinom(100000, 2 * treat.SCL.gaussian.25 + 1, expit(eval(parse(text = base.test)) + treat.SCL.gaussian.25 * (1 + x.test[, 2] + x.test[, 3]^2)))-rbeta(100000, 2, 0.5), 0.25)
  betav.SCL.gaussian.50 <- quantile(rnbinom(100000, 2 * treat.SCL.gaussian.50 + 1, expit(eval(parse(text = base.test)) + treat.SCL.gaussian.50 * (1 + x.test[, 2] + x.test[, 3]^2)))-rbeta(100000, 2, 0.5), 0.50)

  betav.SCL.linear.25.beta <- quantile(rnbinom(100000, 2 * treat.SCL.linear.25.beta + 1, expit(eval(parse(text = base.test)) + treat.SCL.linear.25.beta * (1 + x.test[, 2] + x.test[, 3]^2)))-rbeta(100000, 2, 0.5), 0.25)
  betav.SCL.linear.50.beta <- quantile(rnbinom(100000, 2 * treat.SCL.linear.50.beta + 1, expit(eval(parse(text = base.test)) + treat.SCL.linear.50.beta * (1 + x.test[, 2] + x.test[, 3]^2)))-rbeta(100000, 2, 0.5), 0.50)

  betav.SCL.gaussian.25.beta <- quantile(rnbinom(100000, 2 * treat.SCL.gaussian.25.beta + 1, expit(eval(parse(text = base.test)) + treat.SCL.gaussian.25.beta * (1 + x.test[, 2] + x.test[, 3]^2)))-rbeta(100000, 2, 0.5), 0.25)
  betav.SCL.gaussian.50.beta <- quantile(rnbinom(100000, 2 * treat.SCL.gaussian.50.beta + 1, expit(eval(parse(text = base.test)) + treat.SCL.gaussian.50.beta * (1 + x.test[, 2] + x.test[, 3]^2)))-rbeta(100000, 2, 0.5), 0.50)
  
  return(list(
    SCL.linear = c(v.SCL.linear.25, betav.SCL.linear.25, mr.SCL.linear.25, v.SCL.linear.50, betav.SCL.linear.50, mr.SCL.linear.50),
    SCL.gaussian = c(v.SCL.gaussian.25, betav.SCL.gaussian.25, mr.SCL.gaussian.25, v.SCL.gaussian.50, betav.SCL.gaussian.50, mr.SCL.gaussian.50),
    SCL.linear.beta = c(v.SCL.linear.25.beta, betav.SCL.linear.25.beta, mr.SCL.linear.25.beta, v.SCL.linear.50.beta, betav.SCL.linear.50.beta, mr.SCL.linear.50.beta),
    SCL.gaussian.beta = c(v.SCL.gaussian.25.beta, betav.SCL.gaussian.25.beta, mr.SCL.gaussian.25.beta, v.SCL.gaussian.50.beta, betav.SCL.gaussian.50.beta, mr.SCL.gaussian.50.beta)
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
})
clusterExport(cl, varlist = list("SCL.count", "SCL.count.beta", "smooth.quantile", "expit"), envir = .GlobalEnv)

# parallel computing
system.time({
  res <- parLapply(cl, 1:times, simulation, n = n, d=d, base = case3.base, base.test = case3.base.test, x.test = x, opt.q.25 = case3.opt.q.25, opt.q.50 = case3.opt.q.50)
})
stopCluster(cl)

# extract results
# value funtion and misclassification rates
case3.SCL.linear <- array(0, dim = c(times, 6))
case3.SCL.gaussian <- array(0, dim = c(times, 6))

case3.SCL.linear.beta <- array(0, dim = c(times, 6))
case3.SCL.gaussian.beta <- array(0, dim = c(times, 6))

for (j in 1:times) {
  case3.SCL.linear[j, ] <- res[[j]]$SCL.linear
  case3.SCL.gaussian[j, ] <- res[[j]]$SCL.gaussian

  case3.SCL.linear.beta[j, ] <- res[[j]]$SCL.linear.beta
  case3.SCL.gaussian.beta[j, ] <- res[[j]]$SCL.gaussian.beta
}