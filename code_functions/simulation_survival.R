# import necessary libraries
library(parallel)
library(WeightSVM)
library(survival)

# source necessary functions
source("./code_functions/function_survival.R")

### The main simulation function that generates data, fits different methods, and evaluates their performance.
### The input l denotes the index of the simulation run, the inputs-n, d, C0, base-are used to generate data, and the inputs-base.test, x.test, eps.test, opt.q.25, opt.q.50, q.25, q.50-are used to evaluate different methods. 
###  The outputs include all metrics for competing methods, including the value, the estimated value, the bias of the estimated optimal value, misclassification rates, computation time, and quantile optimal treatment on a grid for plotting purposes
simulation <- function(l, n, d, C0, base, base.test, x.test, eps.test, opt.q.25, opt.q.50, q.25, q.50) {
  # set the sample size and the dimension of the coefficients
  n <- n
  d <- d

  # generate data
  x <- cbind(1, array(rnorm(n * d), c(n, d)))
  ps <- rep(0.5,n) 
  a <- sapply(ps, rbinom, n = 1, size = 1)
  logy <- eval(parse(text = base)) +  a * (x[, 2] + x[, 3]) + (a+1) * rnorm(n)
  y <- exp(logy/4)
  logc <- C0 + 4*a*(x[, 2]+0.625*x[,3]) + rnorm(n)
  c <- exp(logc/4)
  delta <- as.numeric(y <= c)
  tildey <- pmin(y, c)
  data <- list(tildey = tildey, delta = delta, a = a, x = x)

  ### estimation of nuisance functions
  # correct specification of the survival function of the potential outcomes
  data1 <- data.frame(tildey = tildey[a==1], x = x[a==1, -1], delta = delta[a==1])
  data0 <- data.frame(tildey = tildey[a==0], x = x[a==0, -1], delta = delta[a==0])
  f.y.x1.correct <- survreg(Surv(tildey, delta) ~  x.1 + x.2,
              data = data1,
              dist = "lognormal")
  class(f.y.x1.correct) <- "lognorm"
  f.y.x0.correct <- survreg(Surv(tildey, delta) ~  x.1 + x.2,
              data = data0,
              dist = "lognormal")
  class(f.y.x0.correct) <- "lognorm"

  # incorrect specification of the survival function of the potential outcomes
  f.y.x1.incorrect <- survreg(Surv(tildey, delta) ~  x.1 ,
              data = data1,
              dist = "lognormal")
  class(f.y.x1.incorrect) <- "lognorm"
  f.y.x0.incorrect <- survreg(Surv(tildey, delta) ~  x.1 ,
              data = data0,
              dist = "lognormal")
  class(f.y.x0.incorrect) <- "lognorm"

  # correct specification of the survival function of the censoring time 
  data.c <- data.frame(tildey = tildey, x = x[,-1], a = a, delta = delta)
  f.c.ax.correct <- survreg(Surv(tildey, 1-delta) ~ a * x.1 + a * x.2,
              data = data.c,
              dist = "lognormal")
  class(f.c.ax.correct) <- "lognorm"

  # incorrect specification of the survival function of the censoring time
  f.c.ax.incorrect <- survreg(Surv(tildey, 1-delta) ~ a * x.1 + a * I(x.2^3),
              data = data.c,
              dist = "lognormal")
  class(f.c.ax.incorrect) <- "lognorm"

  ### Run competeting methods for estimating quantile OTRs
  # Our successive classification learning (SCL)
  kappa <- sd(tildey) / sqrt(n) / 6
  epsilon <- 1 / sqrt(n) / 2
  
  # correct specification of all nuisance functions
  res.SCL.linear.25.correctcorrect <- SCL.survival(data, 0.25, kappa, epsilon, ps, f.y.x1.correct, f.y.x0.correct, f.c.ax.correct, kernel="linear")
  res.SCL.gaussian.25.correctcorrect <- SCL.survival(data, 0.25, kappa, epsilon, ps, f.y.x1.correct, f.y.x0.correct, f.c.ax.correct, kernel="gaussian")
  res.SCL.linear.50.correctcorrect <- SCL.survival(data, 0.50, kappa, epsilon, ps, f.y.x1.correct, f.y.x0.correct, f.c.ax.correct, kernel="linear")
  res.SCL.gaussian.50.correctcorrect <- SCL.survival(data, 0.50, kappa, epsilon, ps, f.y.x1.correct, f.y.x0.correct, f.c.ax.correct, kernel="gaussian")

  # incorrect specification of censoring model and correct specification of the survival model of the potential outcomes
  res.SCL.linear.25.incorrectcorrect <- SCL.survival(data, 0.25, kappa, epsilon, ps, f.y.x1.correct, f.y.x0.correct, f.c.ax.incorrect, kernel="linear")
  res.SCL.gaussian.25.incorrectcorrect <- SCL.survival(data, 0.25, kappa, epsilon, ps, f.y.x1.correct, f.y.x0.correct, f.c.ax.incorrect, kernel="gaussian")
  res.SCL.linear.50.incorrectcorrect <- SCL.survival(data, 0.50, kappa, epsilon, ps, f.y.x1.correct, f.y.x0.correct, f.c.ax.incorrect, kernel="linear")
  res.SCL.gaussian.50.incorrectcorrect <- SCL.survival(data, 0.50, kappa, epsilon, ps, f.y.x1.correct, f.y.x0.correct, f.c.ax.incorrect, kernel="gaussian")

  # correct specification of the censoring model and incorrect specification of the survival model of the potential outcomes 
  res.SCL.linear.25.correctincorrect <- SCL.survival(data, 0.25, kappa, epsilon, ps, f.y.x1.incorrect, f.y.x0.incorrect, f.c.ax.correct, kernel="linear")
  res.SCL.gaussian.25.correctincorrect <- SCL.survival(data, 0.25, kappa, epsilon, ps, f.y.x1.incorrect, f.y.x0.incorrect, f.c.ax.correct, kernel="gaussian")
  res.SCL.linear.50.correctincorrect <- SCL.survival(data, 0.50, kappa, epsilon, ps, f.y.x1.incorrect, f.y.x0.incorrect, f.c.ax.correct, kernel="linear")
  res.SCL.gaussian.50.correctincorrect <- SCL.survival(data, 0.50, kappa, epsilon, ps, f.y.x1.incorrect, f.y.x0.incorrect, f.c.ax.correct, kernel="gaussian")

  # incorrect specification of all nuisance functions
  res.SCL.linear.25.incorrectincorrect <- SCL.survival(data, 0.25, kappa, epsilon, ps, f.y.x1.incorrect, f.y.x0.incorrect, f.c.ax.incorrect, kernel="linear")
  res.SCL.gaussian.25.incorrectincorrect <- SCL.survival(data, 0.25, kappa, epsilon, ps, f.y.x1.incorrect, f.y.x0.incorrect, f.c.ax.incorrect, kernel="gaussian")
  res.SCL.linear.50.incorrectincorrect <- SCL.survival(data, 0.50, kappa, epsilon, ps, f.y.x1.incorrect, f.y.x0.incorrect, f.c.ax.incorrect, kernel="linear")
  res.SCL.gaussian.50.incorrectincorrect <- SCL.survival(data, 0.50, kappa, epsilon, ps, f.y.x1.incorrect, f.y.x0.incorrect, f.c.ax.incorrect, kernel="gaussian")

  # extract the estimated quantile treatments for competeting methods
  # the estimated quantile treatments for x.test
  # all nuisance functions are correctly specified
  xtest.treat.SCL.linear.25.correctcorrect <- as.numeric(predict(res.SCL.linear.25.correctcorrect$opt.reg, x.test[, -1])) - 1
  xtest.treat.SCL.gaussian.25.correctcorrect <- as.numeric(predict(res.SCL.gaussian.25.correctcorrect$opt.reg, x.test[, -1])) - 1
  xtest.treat.SCL.linear.50.correctcorrect <- as.numeric(predict(res.SCL.linear.50.correctcorrect$opt.reg, x.test[, -1])) - 1
  xtest.treat.SCL.gaussian.50.correctcorrect <- as.numeric(predict(res.SCL.gaussian.50.correctcorrect$opt.reg, x.test[, -1])) - 1

  # incorrect specification of censoring model and correct specification of the survival model of the potential outcomes
  xtest.treat.SCL.linear.25.incorrectcorrect <- as.numeric(predict(res.SCL.linear.25.incorrectcorrect$opt.reg, x.test[, -1])) - 1
  xtest.treat.SCL.gaussian.25.incorrectcorrect <- as.numeric(predict(res.SCL.gaussian.25.incorrectcorrect$opt.reg, x.test[, -1])) - 1
  xtest.treat.SCL.linear.50.incorrectcorrect <- as.numeric(predict(res.SCL.linear.50.incorrectcorrect$opt.reg, x.test[, -1])) - 1
  xtest.treat.SCL.gaussian.50.incorrectcorrect <- as.numeric(predict(res.SCL.gaussian.50.incorrectcorrect$opt.reg, x.test[, -1])) - 1

  # correct specification of the censoring model and incorrect specification of the survival model of the potential outcomes 
  xtest.treat.SCL.linear.25.correctincorrect <- as.numeric(predict(res.SCL.linear.25.correctincorrect$opt.reg, x.test[, -1])) - 1
  xtest.treat.SCL.gaussian.25.correctincorrect <- as.numeric(predict(res.SCL.gaussian.25.correctincorrect$opt.reg, x.test[, -1])) - 1
  xtest.treat.SCL.linear.50.correctincorrect <- as.numeric(predict(res.SCL.linear.50.correctincorrect$opt.reg, x.test[, -1])) - 1
  xtest.treat.SCL.gaussian.50.correctincorrect <- as.numeric(predict(res.SCL.gaussian.50.correctincorrect$opt.reg, x.test[, -1])) - 1

  # incorrect specification of all nuisance functions
  xtest.treat.SCL.linear.25.incorrectincorrect <- as.numeric(predict(res.SCL.linear.25.incorrectincorrect$opt.reg, x.test[, -1])) - 1
  xtest.treat.SCL.gaussian.25.incorrectincorrect <- as.numeric(predict(res.SCL.gaussian.25.incorrectincorrect$opt.reg, x.test[, -1])) - 1
  xtest.treat.SCL.linear.50.incorrectincorrect <- as.numeric(predict(res.SCL.linear.50.incorrectincorrect$opt.reg, x.test[, -1])) - 1
  xtest.treat.SCL.gaussian.50.incorrectincorrect <- as.numeric(predict(res.SCL.gaussian.50.incorrectincorrect$opt.reg, x.test[, -1])) - 1

  ## the bias of the estimated optimal quantile
  # q.SCL.linear.25 <- res.SCL.linear.25$q-q.25
  # q.SCL.gaussian.25 <- res.SCL.gaussian.25$q-q.25
  # q.SCL.linear.50 <- res.SCL.linear.50$q-q.50
  # q.SCL.gaussian.50 <- res.SCL.gaussian.50$q-q.50

  # misspecification rates
  mr.SCL.linear.25.correctcorrect <- 1 - sum(xtest.treat.SCL.linear.25.correctcorrect == opt.q.25) / dim(x.test)[1]
  mr.SCL.gaussian.25.correctcorrect <- 1 - sum(xtest.treat.SCL.gaussian.25.correctcorrect == opt.q.25) / dim(x.test)[1]
  mr.SCL.linear.50.correctcorrect <- 1 - sum(xtest.treat.SCL.linear.50.correctcorrect == opt.q.50) / dim(x.test)[1]
  mr.SCL.gaussian.50.correctcorrect <- 1 - sum(xtest.treat.SCL.gaussian.50.correctcorrect == opt.q.50) / dim(x.test)[1]

  mr.SCL.linear.25.incorrectcorrect <- 1 - sum(xtest.treat.SCL.linear.25.incorrectcorrect == opt.q.25) / dim(x.test)[1]
  mr.SCL.gaussian.25.incorrectcorrect <- 1 - sum(xtest.treat.SCL.gaussian.25.incorrectcorrect == opt.q.25) / dim(x.test)[1]
  mr.SCL.linear.50.incorrectcorrect <- 1 - sum(xtest.treat.SCL.linear.50.incorrectcorrect == opt.q.50) / dim(x.test)[1]
  mr.SCL.gaussian.50.incorrectcorrect <- 1 - sum(xtest.treat.SCL.gaussian.50.incorrectcorrect == opt.q.50) / dim(x.test)[1]

  mr.SCL.linear.25.correctincorrect <- 1 - sum(xtest.treat.SCL.linear.25.correctincorrect == opt.q.25) / dim(x.test)[1]
  mr.SCL.gaussian.25.correctincorrect <- 1 - sum(xtest.treat.SCL.gaussian.25.correctincorrect == opt.q.25) / dim(x.test)[1]
  mr.SCL.linear.50.correctincorrect <- 1 - sum(xtest.treat.SCL.linear.50.correctincorrect == opt.q.50) / dim(x.test)[1]
  mr.SCL.gaussian.50.correctincorrect <- 1 - sum(xtest.treat.SCL.gaussian.50.correctincorrect == opt.q.50) / dim(x.test)[1]

  mr.SCL.linear.25.incorrectincorrect <- 1 - sum(xtest.treat.SCL.linear.25.incorrectincorrect == opt.q.25) / dim(x.test)[1]
  mr.SCL.gaussian.25.incorrectincorrect <- 1 - sum(xtest.treat.SCL.gaussian.25.incorrectincorrect == opt.q.25) / dim(x.test)[1]
  mr.SCL.linear.50.incorrectincorrect <- 1 - sum(xtest.treat.SCL.linear.50.incorrectincorrect == opt.q.50) / dim(x.test)[1]
  mr.SCL.gaussian.50.incorrectincorrect <- 1 - sum(xtest.treat.SCL.gaussian.50.incorrectincorrect == opt.q.50) / dim(x.test)[1]

  # true value function
  v.SCL.linear.25.correctcorrect <-quantile(exp((eval(parse(text = base.test)) + xtest.treat.SCL.linear.25.correctcorrect  * (x.test[, 2] + x.test[, 3]) + ( xtest.treat.SCL.linear.25.correctcorrect + 1) * eps.test)/4), 0.25)
  v.SCL.linear.50.correctcorrect <- quantile(exp((eval(parse(text = base.test)) + xtest.treat.SCL.linear.50.correctcorrect  * (x.test[, 2] + x.test[, 3]) + (xtest.treat.SCL.linear.50.correctcorrect + 1) * eps.test)/4), 0.50)
  v.SCL.gaussian.25.correctcorrect <- quantile(exp((eval(parse(text = base.test)) + xtest.treat.SCL.gaussian.25.correctcorrect  * (x.test[, 2] + x.test[, 3]) + ( xtest.treat.SCL.gaussian.25.correctcorrect + 1) * eps.test)/4), 0.25)
  v.SCL.gaussian.50.correctcorrect <- quantile(exp((eval(parse(text = base.test)) + xtest.treat.SCL.gaussian.50.correctcorrect  * (x.test[, 2] + x.test[, 3]) + ( xtest.treat.SCL.gaussian.50.correctcorrect + 1) * eps.test)/4), 0.50)

  v.SCL.gaussian.25.incorrectcorrect <- quantile(exp((eval(parse(text = base.test)) + xtest.treat.SCL.gaussian.25.incorrectcorrect  * (x.test[, 2] + x.test[, 3]) + ( xtest.treat.SCL.gaussian.25.incorrectcorrect + 1) * eps.test)/4), 0.25)
  v.SCL.gaussian.50.incorrectcorrect <- quantile(exp((eval(parse(text = base.test)) + xtest.treat.SCL.gaussian.50.incorrectcorrect  * (x.test[, 2] + x.test[, 3]) + ( xtest.treat.SCL.gaussian.50.incorrectcorrect + 1) * eps.test)/4), 0.50)
  v.SCL.linear.25.incorrectcorrect <- quantile(exp((eval(parse(text = base.test)) + xtest.treat.SCL.linear.25.incorrectcorrect  * (x.test[, 2] + x.test[, 3]) + ( xtest.treat.SCL.linear.25.incorrectcorrect + 1) * eps.test)/4), 0.25)
  v.SCL.linear.50.incorrectcorrect <- quantile(exp((eval(parse(text = base.test)) + xtest.treat.SCL.linear.50.incorrectcorrect  * (x.test[, 2] + x.test[, 3]) + ( xtest.treat.SCL.linear.50.incorrectcorrect + 1) * eps.test)/4), 0.50)

  v.SCL.linear.25.correctincorrect <- quantile(exp((eval(parse(text = base.test)) + xtest.treat.SCL.linear.25.correctincorrect  * (x.test[, 2] + x.test[, 3]) + ( xtest.treat.SCL.linear.25.correctincorrect + 1) * eps.test)/4), 0.25)
  v.SCL.linear.50.correctincorrect <- quantile(exp((eval(parse(text = base.test)) + xtest.treat.SCL.linear.50.correctincorrect  * (x.test[, 2] + x.test[, 3]) + ( xtest.treat.SCL.linear.50.correctincorrect + 1) * eps.test)/4), 0.50)
  v.SCL.gaussian.25.correctincorrect <- quantile(exp((eval(parse(text = base.test)) + xtest.treat.SCL.gaussian.25.correctincorrect  * (x.test[, 2] + x.test[, 3]) + ( xtest.treat.SCL.gaussian.25.correctincorrect + 1) * eps.test)/4), 0.25)
  v.SCL.gaussian.50.correctincorrect <- quantile(exp((eval(parse(text = base.test)) + xtest.treat.SCL.gaussian.50.correctincorrect  * (x.test[, 2] + x.test[, 3]) + ( xtest.treat.SCL.gaussian.50.correctincorrect + 1) * eps.test)/4), 0.50)

  v.SCL.linear.25.incorrectincorrect <- quantile(exp((eval(parse(text = base.test)) + xtest.treat.SCL.linear.25.incorrectincorrect  * (x.test[, 2] + x.test[, 3]) + ( xtest.treat.SCL.linear.25.incorrectincorrect + 1) * eps.test)/4), 0.25)
  v.SCL.linear.50.incorrectincorrect <- quantile(exp((eval(parse(text = base.test)) + xtest.treat.SCL.linear.50.incorrectincorrect  * (x.test[, 2] + x.test[, 3]) + ( xtest.treat.SCL.linear.50.incorrectincorrect + 1) * eps.test)/4), 0.50)
  v.SCL.gaussian.25.incorrectincorrect <- quantile(exp((eval(parse(text = base.test)) + xtest.treat.SCL.gaussian.25.incorrectincorrect  * (x.test[, 2] + x.test[, 3]) + ( xtest.treat.SCL.gaussian.25.incorrectincorrect + 1) * eps.test)/4), 0.25)
  v.SCL.gaussian.50.incorrectincorrect <- quantile(exp((eval(parse(text = base.test)) + xtest.treat.SCL.gaussian.50.incorrectincorrect  * (x.test[, 2] + x.test[, 3]) + ( xtest.treat.SCL.gaussian.50.incorrectincorrect + 1) * eps.test)/4), 0.50)

  return(list(
    SCL.linear.correctcorrect = c(v.SCL.linear.25.correctcorrect, mr.SCL.linear.25.correctcorrect, v.SCL.linear.50.correctcorrect, mr.SCL.linear.50.correctcorrect),
    SCL.gaussian.correctcorrect = c(v.SCL.gaussian.25.correctcorrect, mr.SCL.gaussian.25.correctcorrect, v.SCL.gaussian.50.correctcorrect, mr.SCL.gaussian.50.correctcorrect),
    SCL.linear.incorrectcorrect = c(v.SCL.linear.25.incorrectcorrect, mr.SCL.linear.25.incorrectcorrect, v.SCL.linear.50.incorrectcorrect, mr.SCL.linear.50.incorrectcorrect),
    SCL.gaussian.incorrectcorrect = c(v.SCL.gaussian.25.incorrectcorrect, mr.SCL.gaussian.25.incorrectcorrect, v.SCL.gaussian.50.incorrectcorrect, mr.SCL.gaussian.50.incorrectcorrect),
    SCL.linear.correctincorrect = c(v.SCL.linear.25.correctincorrect, mr.SCL.linear.25.correctincorrect, v.SCL.linear.50.correctincorrect, mr.SCL.linear.50.correctincorrect),
    SCL.gaussian.correctincorrect = c(v.SCL.gaussian.25.correctincorrect, mr.SCL.gaussian.25.correctincorrect, v.SCL.gaussian.50.correctincorrect, mr.SCL.gaussian.50.correctincorrect),
    SCL.linear.incorrectincorrect = c(v.SCL.linear.25.incorrectincorrect, mr.SCL.linear.25.incorrectincorrect, v.SCL.linear.50.incorrectincorrect, mr.SCL.linear.50.incorrectincorrect),
    SCL.gaussian.incorrectincorrect = c(v.SCL.gaussian.25.incorrectincorrect, mr.SCL.gaussian.25.incorrectincorrect, v.SCL.gaussian.50.incorrectincorrect, mr.SCL.gaussian.50.incorrectincorrect)
  ))
}

# Specifiy the data generating process. If you want to test other data generating processes, please change the base, n, err, err.test accordingly.
# Note that the first column of x is the intercept 1.
set.seed(0)
case1.base <- "1.25-1*x[,2]-0.5*x[,3]"
case1.base.test <-"1.25-1*x.test[,2]-0.5*x.test[,3]"
case1.err <- "rnorm(n)"
case1.err.test <- "rnorm(m)"

# produce testing set
m <- 100000 # sample size
x <- cbind(1, array(rnorm(m * d), c(m, d)))
ps <- rep(0.5, m)
a <- sapply(ps, rbinom, n = 1, size = 1)
eps <- eval(parse(text = case1.err.test))
logy <- eval(parse(text = case1.base)) +  a * (x[, 2] + x[, 3]) + (a+1) * eps
y <- exp(logy/4)
logc <- C0 + 4*a*(x[, 2]+0.625*x[,3]) + rnorm(m)
c <- exp(logc/4)
tildey <- pmin(y, c)
delta <- as.numeric(y <= c)

# calculate the true optimal quantile 
case1.q.25 <- uniroot(opt.q.case1, c(min(tildey), 100), tau = 0.25)$root
case1.q.50 <- uniroot(opt.q.case1, c(min(tildey), 100), tau = 0.50)$root

# calculate the true optimal treatment regime for testing set
case1.opt.q.25 <- opt.reg.case1(case1.q.25)
case1.opt.q.50 <- opt.reg.case1(case1.q.50)

# parallel setup
cl <- makeCluster(getOption("cl.cores", clnum))
clusterSetRNGStream(cl, 999)
clusterEvalQ(cl, {
  library(WeightSVM)
  library(survival)
})
clusterExport(cl, varlist = list("SCL.survival", "predict.lognorm"), envir = .GlobalEnv)

# parallel computing
system.time({
  res <- parLapply(cl, 1:times, simulation, n = n, d=d, C0 = C0, base = case1.base, base.test = case1.base.test, x.test = x, eps.test = eps, opt.q.25 = case1.opt.q.25, opt.q.50 = case1.opt.q.50, q.25 = case1.q.25, q.50 = case1.q.50)
})
stopCluster(cl)

# extract results
# value funtion and misclassification rates
case1.SCL.linear.correctcorrect <- array(0, dim = c(times, 4))
case1.SCL.gaussian.correctcorrect <- array(0, dim = c(times, 4))
case1.SCL.linear.incorrectcorrect <- array(0, dim = c(times, 4))
case1.SCL.gaussian.incorrectcorrect <- array(0, dim = c(times, 4))
case1.SCL.linear.correctincorrect <- array(0, dim = c(times, 4))
case1.SCL.gaussian.correctincorrect <- array(0, dim = c(times, 4))
case1.SCL.linear.incorrectincorrect <- array(0, dim = c(times, 4))
case1.SCL.gaussian.incorrectincorrect <- array(0, dim = c(times, 4))

for (j in 1:times) {
  case1.SCL.linear.correctcorrect[j, ] <- res[[j]]$SCL.linear.correctcorrect
  case1.SCL.gaussian.correctcorrect[j, ] <- res[[j]]$SCL.gaussian.correctcorrect
  case1.SCL.linear.incorrectcorrect[j, ] <- res[[j]]$SCL.linear.incorrectcorrect
  case1.SCL.gaussian.incorrectcorrect[j, ] <- res[[j]]$SCL.gaussian.incorrectcorrect
  case1.SCL.linear.correctincorrect[j, ] <- res[[j]]$SCL.linear.correctincorrect
  case1.SCL.gaussian.correctincorrect[j, ] <- res[[j]]$SCL.gaussian.correctincorrect
  case1.SCL.linear.incorrectincorrect[j, ] <- res[[j]]$SCL.linear.incorrectincorrect
  case1.SCL.gaussian.incorrectincorrect[j, ] <- res[[j]]$SCL.gaussian.incorrectincorrect
}
