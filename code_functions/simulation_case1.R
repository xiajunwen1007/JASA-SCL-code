# import necessary libraries
library(MASS)
library(parallel)
library(glmnet)
library(WeightSVM)
library(quantoptr)
library(quantreg)

# source necessary functions
source("./code_functions/function_main.R")

### The main simulation function that generates data, fits different methods, and evaluates their performance.
### The input l denotes the index of the simulation run, the inputs-n, d, err, base-are used to generate data, and the inputs-base.test, x.test, eps.test, opt.q.25, opt.q.50, q.25, q.50-are used to evaluate different methods. 
###  The outputs include all metrics for competing methods, including the value, the estimated value, the bias of the estimated optimal value, misclassification rates, computation time, and quantile optimal treatment on a grid for plotting purposes
simulation <- function(l, n, d, err, base, base.test, x.test, eps.test, opt.q.25, opt.q.50, q.25, q.50) {
  # set the sample size and the dimension of the coefficients
  n <- n
  d <- d

  # generate data
  x <- cbind(1, array(rnorm(n * d), c(n, d)))
  ps <- exp(-0.5 - 0.5 * (x[, 2] + x[, 3])) / (1 + exp(-0.5 - 0.5 * (x[, 2] + x[, 3])))
  a <- sapply(ps, rbinom, n = 1, size = 1)
  eps <- eval(parse(text = err))
  y <- eval(parse(text = base)) +  a * (x[, 2] + x[, 3]) + (2 * a + 1) * eps
  data <- list(y = y, x = x, a = a)

  time.nuisance.ps.0 <- Sys.time()
  # estimation of the propensity score model
  if(d<=5){
    f.a.x <- glmnet(x[, -1], a, family = binomial(link = "logit"), lambda = 0)
    mean.a <- predict(f.a.x, x[, -1], type = "response")
  }else{
    f.a.x <- cv.glmnet(x[, -1], a, family = binomial(link = "logit"))
    mean.a <- predict(f.a.x, x[, -1], type = "response")
  }
  time.nuisance.ps.1 <- Sys.time()

  time.nuisance.sm.0 <- Sys.time()
  # estimation of the conditional survival model
  dataframe <- as.data.frame(list(y=y,x=x[,-1],a=a))
  formula <-"y ~ "
  for(i in 1:d){
    formula = paste0(formula, paste("x.", i, sep = ""))
    formula = paste0(formula, " + ")
    formula = paste0(formula, paste("I(a*x.", i, ")", sep = ""))
    formula = paste0(formula, " + ")
  }
  formula <- as.formula(paste0(formula, "a"))
  xa <- model.matrix(formula, dataframe)[, -1]

  if(d<=5){
    fit.f.y.ax <- glmnet(xa, y, lambda=0)
    f.y.ax <- list(fit=fit.f.y.ax, formula=formula)
    class(f.y.ax) <- "glmnetformula"

    r <- y - predict(fit.f.y.ax, xa)
    fit.f.r.ax <- glmnet(xa, r^2, family = Gamma(link = "log"), lambda=0)
    f.r.ax <- list(fit=fit.f.r.ax, formula=formula)
    class(f.r.ax) <- "glmnetformula"
  }else{
    fit.f.y.ax <- cv.glmnet(xa, y)
    f.y.ax <- list(fit=fit.f.y.ax, formula=formula)
    class(f.y.ax) <- "glmnetformula"

    r <- y - predict(fit.f.y.ax, xa)
    fit.f.r.ax <- cv.glmnet(xa, r^2, family = Gamma(link = "log"))
    f.r.ax <- list(fit=fit.f.r.ax, formula=formula)
    class(f.r.ax) <- "glmnetformula"
  }
  

  time.nuisance.sm.1 <- Sys.time()

  ### Run competeting methods for estimating quantile OTRs
  # Our successive classification learning (SCL)
  kappa <- sd(y) / sqrt(n) / 6
  epsilon <- 1 / sqrt(n) / 2

  time_0 <- Sys.time()
  res.SCL.linear.25 <- SCL.continuous(data, 0.25, kappa, epsilon, mean.a, f.y.ax, f.r.ax, kernel="linear")
  time_1 <- Sys.time()
  res.SCL.gaussian.25 <- SCL.continuous(data, 0.25, kappa, epsilon, mean.a, f.y.ax, f.r.ax, kernel="gaussian")
  time_2 <- Sys.time()
  res.SCL.linear.50 <- SCL.continuous(data, 0.50, kappa, epsilon, mean.a, f.y.ax, f.r.ax, kernel="linear")
  time_3 <- Sys.time()
  res.SCL.gaussian.50 <- SCL.continuous(data, 0.50, kappa, epsilon, mean.a, f.y.ax, f.r.ax, kernel="gaussian")
  time_4 <- Sys.time()

  # QIQ learning
  res.QIQ.25 <- QIQ.learning.continuous(data, 0.25, f.y.ax, f.r.ax)
  time_5 <- Sys.time()
  res.QIQ.50 <- QIQ.learning.continuous(data, 0.50, f.y.ax, f.r.ax)
  time_6 <- Sys.time()

  # Wang's doubly robust method
  xformula <- ""
  for(j in 1:(d)){
    xformula <- paste(xformula,paste("x.", j, sep = ""),seq="")
    if (j < d) {
      xformula <- paste(xformula,"+",sep="")
    }
  }
  data.wang <- data.frame(y = y, x = x[, -1], a = a)
  regimeClass <- as.formula(paste("a ~", xformula))
  moCondQuant_0 <- as.formula(y~x.1+x.2)
  moCondQuant_1 <- as.formula(y~x.1+x.2)
  moPropen <- as.formula(a~x.1+x.2)

  time_7 <- Sys.time()
  res.wang.25 <- DR_Qopt(
    data = data.wang, regimeClass = regimeClass, tau = 0.25,
    moPropen = moPropen,
    moCondQuant_0 = moCondQuant_0,
    moCondQuant_1 = moCondQuant_1,
    pop.size = 500, it.num = 2
  )
  time_8 <- Sys.time()
  res.wang.50 <- DR_Qopt(
    data = data.wang, regimeClass = regimeClass, tau = 0.5,
    moPropen = moPropen,
    moCondQuant_0 = moCondQuant_0,
    moCondQuant_1 = moCondQuant_1,
    pop.size = 500, it.num = 2
  )
  time_9 <- Sys.time()

  # SCL with hard regime. Note that this method is only evaluated by the value function and the misspecification rate
  res.SCL.linear.25.hard <- SCL.continuous.hard(data, 0.25, kappa, epsilon, mean.a, f.y.ax, f.r.ax, kernel = "linear")
  res.SCL.gaussian.25.hard <- SCL.continuous.hard(data, 0.25, kappa, epsilon, mean.a, f.y.ax, f.r.ax, kernel = "gaussian")
  res.SCL.linear.50.hard <- SCL.continuous.hard(data, 0.50, kappa, epsilon, mean.a, f.y.ax, f.r.ax, kernel = "linear")
  res.SCL.gaussian.50.hard <- SCL.continuous.hard(data, 0.50, kappa, epsilon, mean.a, f.y.ax, f.r.ax, kernel = "gaussian")

  # extract the estimated quantile treatments for competeting methods
  # the estimated quantile treatments for x.test and x
  xtest.treat.SCL.linear.25 <- as.numeric(predict(res.SCL.linear.25$opt.reg, x.test[, -1])) - 1
  xtreat.SCL.linear.25 <- as.numeric(predict(res.SCL.linear.25$opt.reg, x[, -1]))-1

  xtest.treat.SCL.gaussian.25 <- as.numeric(predict(res.SCL.gaussian.25$opt.reg, x.test[, -1])) - 1
  xtreat.SCL.gaussian.25 <- as.numeric(predict(res.SCL.gaussian.25$opt.reg, x[, -1]))-1

  xtest.treat.SCL.linear.50 <- as.numeric(predict(res.SCL.linear.50$opt.reg, x.test[, -1])) - 1
  xtreat.SCL.linear.50 <- as.numeric(predict(res.SCL.linear.50$opt.reg, x[, -1]))-1

  xtest.treat.SCL.gaussian.50 <- as.numeric(predict(res.SCL.gaussian.50$opt.reg, x.test[, -1])) - 1
  xtreat.SCL.gaussian.50 <- as.numeric(predict(res.SCL.gaussian.50$opt.reg, x[, -1]))-1

  xtest.treat.QIQ.25 <- predict(res.QIQ.25$opt.reg, x.test[,-1])
  xtreat.QIQ.25 <- predict(res.QIQ.25$opt.reg, x[,-1])

  xtest.treat.QIQ.50 <- predict(res.QIQ.50$opt.reg, x.test[,-1])
  xtreat.QIQ.50 <- predict(res.QIQ.50$opt.reg, x[,-1])

  if (sum(res.wang.25$coefficients) == "NaN") {
    xtest.treat.wang.25 <- x.test %*% rep(0, d+1)
    xtreat.wang.25 <- x %*% rep(0, d+1)
  } else {
    xtest.treat.wang.25 <- x.test %*% res.wang.25$coefficients > 0
    xtreat.wang.25 <- x %*% res.wang.25$coefficients > 0
  }

  if (sum(res.wang.50$coefficients) == "NaN") {
    xtest.treat.wang.50 <- x.test %*% rep(0, d+1)
    xtreat.wang.50 <- x %*% rep(0, d+1)
  } else {
    xtest.treat.wang.50 <- x.test %*% res.wang.50$coefficients > 0
    xtreat.wang.50 <- x %*% res.wang.50$coefficients > 0
  }

  xtest.treat.SCL.linear.25.hard <- as.numeric(predict(res.SCL.linear.25.hard$opt.reg, x.test[, -1])) - 1
  xtest.treat.SCL.gaussian.25.hard <- as.numeric(predict(res.SCL.gaussian.25.hard$opt.reg, x.test[, -1])) - 1
  xtest.treat.SCL.linear.50.hard <- as.numeric(predict(res.SCL.linear.50.hard$opt.reg, x.test[, -1])) - 1
  xtest.treat.SCL.gaussian.50.hard <- as.numeric(predict(res.SCL.gaussian.50.hard$opt.reg, x.test[, -1])) - 1

  ###Evaluate competing methods
  # time consuption
  time.nuisance.ps <- difftime(time.nuisance.ps.1, time.nuisance.ps.0, units = "secs")
  time.nuisance.sm <- difftime(time.nuisance.sm.1, time.nuisance.sm.0, units = "secs")

  time.SCL.linear.25 <- difftime(time_1, time_0, units = "secs")+ time.nuisance.ps+time.nuisance.sm 
  time.SCL.gaussian.25 <- difftime(time_2, time_1, units = "secs")+  time.nuisance.ps+time.nuisance.sm 
  time.SCL.linear.50 <- difftime(time_3, time_2, units = "secs")+  time.nuisance.ps+time.nuisance.sm 
  time.SCL.gaussian.50 <- difftime(time_4, time_3, units = "secs")+  time.nuisance.ps+time.nuisance.sm 
  time.QIQ.25 <- difftime(time_5, time_4, units = "secs") + time.nuisance.sm
  time.QIQ.50 <- difftime(time_6, time_5, units = "secs") + time.nuisance.sm
  time.wang.25 <- difftime(time_8, time_7, units = "secs") 
  time.wang.50 <- difftime(time_9, time_8, units = "secs")

  # the bias of the estimated optimal quantile
  q.SCL.linear.25 <- res.SCL.linear.25$q-q.25
  q.SCL.gaussian.25 <- res.SCL.gaussian.25$q-q.25
  q.SCL.linear.50 <- res.SCL.linear.50$q-q.50
  q.SCL.gaussian.50 <- res.SCL.gaussian.50$q-q.50
  q.QIQ.25 <- res.QIQ.25$q-q.25
  q.QIQ.50 <- res.QIQ.50$q-q.50
  q.wang.25 <- res.wang.25$hatQ-q.25
  q.wang.50 <- res.wang.50$hatQ-q.50

  # misspecification rates
  mr.SCL.linear.25 <- 1 - sum(xtest.treat.SCL.linear.25 == opt.q.25) / dim(x.test)[1]
  mr.SCL.gaussian.25 <- 1 - sum(xtest.treat.SCL.gaussian.25 == opt.q.25) / dim(x.test)[1]
  mr.SCL.linear.50 <- 1 - sum(xtest.treat.SCL.linear.50 == opt.q.50) / dim(x.test)[1]
  mr.SCL.gaussian.50 <- 1 - sum(xtest.treat.SCL.gaussian.50 == opt.q.50) / dim(x.test)[1]
  mr.QIQ.25 <- 1 - sum(xtest.treat.QIQ.25 == opt.q.25) / dim(x.test)[1]
  mr.QIQ.50 <- 1 - sum(xtest.treat.QIQ.50 == opt.q.50) / dim(x.test)[1]
  mr.wang.25 <- 1 - sum(xtest.treat.wang.25 == opt.q.25) / dim(x.test)[1]
  mr.wang.50 <- 1 - sum(xtest.treat.wang.50 == opt.q.50) / dim(x.test)[1]

  mr.SCL.linear.25.hard <- 1 - sum(xtest.treat.SCL.linear.25.hard == opt.q.25) / dim(x.test)[1]
  mr.SCL.gaussian.25.hard <- 1 - sum(xtest.treat.SCL.gaussian.25.hard == opt.q.25) / dim(x.test)[1]
  mr.SCL.linear.50.hard <- 1 - sum(xtest.treat.SCL.linear.50.hard == opt.q.50) / dim(x.test)[1]
  mr.SCL.gaussian.50.hard <- 1 - sum(xtest.treat.SCL.gaussian.50.hard == opt.q.50) / dim(x.test)[1]

  # true value function
  v.SCL.linear.25 <- quantile(eval(parse(text = base.test)) + xtest.treat.SCL.linear.25  * (x.test[, 2] + x.test[, 3]) + (2 * xtest.treat.SCL.linear.25 + 1) * eps.test, 0.25)
  v.SCL.linear.50 <- quantile(eval(parse(text = base.test)) + xtest.treat.SCL.linear.50  * (x.test[, 2] + x.test[, 3]) + (2 * xtest.treat.SCL.linear.50 + 1) * eps.test, 0.50)

  v.SCL.gaussian.25 <- quantile(eval(parse(text = base.test)) + xtest.treat.SCL.gaussian.25  * (x.test[, 2] + x.test[, 3]) + (2 * xtest.treat.SCL.gaussian.25 + 1) * eps.test, 0.25)
  v.SCL.gaussian.50 <- quantile(eval(parse(text = base.test)) + xtest.treat.SCL.gaussian.50  * (x.test[, 2] + x.test[, 3]) + (2 * xtest.treat.SCL.gaussian.50 + 1) * eps.test, 0.50)

  v.QIQ.25 <- quantile(eval(parse(text = base.test)) + xtest.treat.QIQ.25  * (x.test[, 2] + x.test[, 3]) + (2 * xtest.treat.QIQ.25 + 1) * eps.test, 0.25)
  v.QIQ.50 <- quantile(eval(parse(text = base.test)) + xtest.treat.QIQ.50  * (x.test[, 2] + x.test[, 3]) + (2 * xtest.treat.QIQ.50 + 1) * eps.test, 0.50)

  v.wang.25 <- quantile(eval(parse(text = base.test)) + xtest.treat.wang.25 * (x.test[, 2] + x.test[, 3]) + (2 * xtest.treat.wang.25 + 1) * eps.test, 0.25)
  v.wang.50 <- quantile(eval(parse(text = base.test)) + xtest.treat.wang.50  * (x.test[, 2] + x.test[, 3]) + (2 * xtest.treat.wang.50 + 1) * eps.test, 0.50)

  v.SCL.linear.25.hard <- quantile(eval(parse(text = base.test)) + xtest.treat.SCL.linear.25.hard * (x.test[, 2] + x.test[, 3]) + (2 * xtest.treat.SCL.linear.25.hard + 1) * eps.test, 0.25)
  v.SCL.linear.50.hard <- quantile(eval(parse(text = base.test)) + xtest.treat.SCL.linear.50.hard * (x.test[, 2] + x.test[, 3]) + (2 * xtest.treat.SCL.linear.50.hard + 1) * eps.test, 0.50)
  v.SCL.gaussian.25.hard <- quantile(eval(parse(text = base.test)) + xtest.treat.SCL.gaussian.25.hard * (x.test[, 2] + x.test[, 3]) + (2 * xtest.treat.SCL.gaussian.25.hard + 1) * eps.test, 0.25)
  v.SCL.gaussian.50.hard <- quantile(eval(parse(text = base.test)) + xtest.treat.SCL.gaussian.50.hard * (x.test[, 2] + x.test[, 3]) + (2 * xtest.treat.SCL.gaussian.50.hard + 1) * eps.test, 0.50)

  ## estimated value function
  hatv.SCL.linear.25 <- rq(y ~ 1, tau = 0.25, weights=as.numeric(xtreat.SCL.linear.25==a)/(a*mean.a+(1-a)*(1-mean.a)))$coefficients
  hatv.SCL.linear.50 <- rq(y ~ 1, tau = 0.50, weights=as.numeric(xtreat.SCL.linear.50==a)/(a*mean.a+(1-a)*(1-mean.a)))$coefficients

  hatv.SCL.gaussian.25 <- rq(y ~ 1, tau = 0.25, weights=as.numeric(xtreat.SCL.gaussian.25==a)/(a*mean.a+(1-a)*(1-mean.a)))$coefficients
  hatv.SCL.gaussian.50 <- rq(y ~ 1, tau = 0.50, weights=as.numeric(xtreat.SCL.gaussian.50==a)/(a*mean.a+(1-a)*(1-mean.a)))$coefficients

  hatv.QIQ.25 <- rq(y ~ 1, tau = 0.25, weights=as.numeric(xtreat.QIQ.25==a)/(a*mean.a+(1-a)*(1-mean.a)))$coefficients
  hatv.QIQ.50 <- rq(y ~ 1, tau = 0.50, weights=as.numeric(xtreat.QIQ.50==a)/(a*mean.a+(1-a)*(1-mean.a)))$coefficients

  hatv.wang.25 <- rq(y ~ 1, tau = 0.25, weights=as.numeric(xtreat.wang.25==a)/(a*mean.a+(1-a)*(1-mean.a)))$coefficients
  hatv.wang.50 <- rq(y ~ 1, tau = 0.50, weights=as.numeric(xtreat.wang.50==a)/(a*mean.a+(1-a)*(1-mean.a)))$coefficients

  true.25 <- as.numeric((q.25 - eval(parse(text = base))) / 1 - ((q.25 - eval(parse(text = base)) - x[, 2] - x[, 3]) / 3) > 0)
  true.50 <- as.numeric((q.50 - eval(parse(text = base))) / 1 - ((q.50 - eval(parse(text = base)) - x[, 2] - x[, 3]) / 3) > 0)

  hatv.true.25 <- rq(y ~ 1, tau = 0.25, weights=as.numeric(true.25==a)/(a*mean.a+(1-a)*(1-mean.a)))$coefficients
  hatv.true.50 <- rq(y ~ 1, tau = 0.50, weights=as.numeric(true.50==a)/(a*mean.a+(1-a)*(1-mean.a)))$coefficients

  # Calculate quantile optimal treatment on a grid for plotting purpose, which is only feasible when d=2.
  if (d==2){
    df <- expand.grid(x1 = seq(-2, 2, length.out = 500),
                    x2 = seq(-2, 2, length.out = 500))
    x.grid <- cbind(1,df,deparse.level = 0)
    names(x.grid)[2] <- "x.1"
    names(x.grid)[3] <- "x.2"

    xgrid.treat.SCL.linear.25 <- as.numeric(predict(res.SCL.linear.25$opt.reg, x.grid[,-1])) - 1
    xgrid.treat.SCL.gaussian.25 <- as.numeric(predict(res.SCL.gaussian.25$opt.reg, x.grid[,-1])) - 1
    xgrid.treat.SCL.linear.50 <- as.numeric(predict(res.SCL.linear.50$opt.reg, x.grid[,-1])) - 1
    xgrid.treat.SCL.gaussian.50 <- as.numeric(predict(res.SCL.gaussian.50$opt.reg, x.grid[,-1])) - 1
    xgrid.treat.QIQ.25 <- as.numeric(predict(res.QIQ.25$opt.reg, unname(as.matrix(x.grid)[,-1])))
    xgrid.treat.QIQ.50 <- as.numeric(predict(res.QIQ.50$opt.reg, unname(as.matrix(x.grid)[,-1])))

    if (sum(res.wang.25$coefficients) == "NaN") {
      xgrid.treat.wang.25 <- as.matrix(x.grid) %*% rep(0, d+1)
    } else {
      xgrid.treat.wang.25 <- as.numeric(as.matrix(x.grid) %*% res.wang.25$coefficients > 0)
    }

    if (sum(res.wang.50$coefficients) == "NaN") {
      xgrid.treat.wang.50 <- as.matrix(x.grid) %*% rep(0, d+1)
    } else {
      xgrid.treat.wang.50 <- as.numeric(as.matrix(x.grid) %*% res.wang.50$coefficients > 0)
    }

    return(list(
      SCL.linear = c(v.SCL.linear.25, mr.SCL.linear.25, v.SCL.linear.50, mr.SCL.linear.50),
      SCL.gaussian = c(v.SCL.gaussian.25, mr.SCL.gaussian.25, v.SCL.gaussian.50, mr.SCL.gaussian.50),
      QIQ = c(v.QIQ.25, mr.QIQ.25, v.QIQ.50, mr.QIQ.50),
      wang = c(v.wang.25, mr.wang.25, v.wang.50, mr.wang.50),
      SCL.linear.hard = c(v.SCL.linear.25.hard, mr.SCL.linear.25.hard, v.SCL.linear.50.hard, mr.SCL.linear.50.hard),
      SCL.gaussian.hard = c(v.SCL.gaussian.25.hard, mr.SCL.gaussian.25.hard, v.SCL.gaussian.50.hard, mr.SCL.gaussian.50.hard),
      time.SCL.linear = c(time.SCL.linear.25, time.SCL.linear.50),
      time.SCL.gaussian = c(time.SCL.gaussian.25, time.SCL.gaussian.50),
      time.QIQ = c(time.QIQ.25, time.QIQ.50),
      time.wang = c(time.wang.25, time.wang.50),
      q.SCL.linear = c(q.SCL.linear.25, q.SCL.linear.50),
      q.SCL.gaussian = c(q.SCL.gaussian.25, q.SCL.gaussian.50),
      q.QIQ = c(q.QIQ.25, q.QIQ.50),
      q.wang = c(q.wang.25, q.wang.50),
      hatv.SCL.linear= c(hatv.SCL.linear.25, hatv.SCL.linear.50),
      hatv.SCL.gaussian= c(hatv.SCL.gaussian.25, hatv.SCL.gaussian.50),
      hatv.QIQ= c(hatv.QIQ.25, hatv.QIQ.50),
      hatv.wang= c(hatv.wang.25, hatv.wang.50),
      hatv.true= c(hatv.true.25, hatv.true.50),
      xgrid.treat.SCL.linear.25 = xgrid.treat.SCL.linear.25,
      xgrid.treat.SCL.linear.50 = xgrid.treat.SCL.linear.50,
      xgrid.treat.SCL.gaussian.25 = xgrid.treat.SCL.gaussian.25,
      xgrid.treat.SCL.gaussian.50 = xgrid.treat.SCL.gaussian.50,
      xgrid.treat.QIQ.25 = xgrid.treat.QIQ.25,
      xgrid.treat.QIQ.50 = xgrid.treat.QIQ.50,
      xgrid.treat.wang.25 = xgrid.treat.wang.25,
      xgrid.treat.wang.50 = xgrid.treat.wang.50
    )) 
  }else{
    return(list(
      SCL.linear = c(v.SCL.linear.25, mr.SCL.linear.25, v.SCL.linear.50, mr.SCL.linear.50),
      SCL.gaussian = c(v.SCL.gaussian.25, mr.SCL.gaussian.25, v.SCL.gaussian.50, mr.SCL.gaussian.50),
      QIQ = c(v.QIQ.25, mr.QIQ.25, v.QIQ.50, mr.QIQ.50),
      wang = c(v.wang.25, mr.wang.25, v.wang.50, mr.wang.50),
      SCL.linear.hard = c(v.SCL.linear.25.hard, mr.SCL.linear.25.hard, v.SCL.linear.50.hard, mr.SCL.linear.50.hard),
      SCL.gaussian.hard = c(v.SCL.gaussian.25.hard, mr.SCL.gaussian.25.hard, v.SCL.gaussian.50.hard, mr.SCL.gaussian.50.hard),
      time.SCL.linear = c(time.SCL.linear.25, time.SCL.linear.50),
      time.SCL.gaussian = c(time.SCL.gaussian.25, time.SCL.gaussian.50),
      time.QIQ = c(time.QIQ.25, time.QIQ.50),
      time.wang = c(time.wang.25, time.wang.50),
      q.SCL.linear = c(q.SCL.linear.25, q.SCL.linear.50),
      q.SCL.gaussian = c(q.SCL.gaussian.25, q.SCL.gaussian.50),
      q.QIQ = c(q.QIQ.25, q.QIQ.50),
      q.wang = c(q.wang.25, q.wang.50),
      hatv.SCL.linear= c(hatv.SCL.linear.25, hatv.SCL.linear.50),
      hatv.SCL.gaussian= c(hatv.SCL.gaussian.25, hatv.SCL.gaussian.50),
      hatv.QIQ= c(hatv.QIQ.25, hatv.QIQ.50),
      hatv.wang= c(hatv.wang.25, hatv.wang.50),
      hatv.true= c(hatv.true.25, hatv.true.50)
    ))
  }
}

# Specifiy the data generating process. If you want to test other data generating processes, please change the base, n, err, err.test accordingly.
# Note that the first column of x is the intercept 1.
set.seed(0)
case1.base <- "2-2*x[,2]+2*x[,3]"
case1.base.test <-"2-2*x.test[,2]+2*x.test[,3]"
case1.err <- "rnorm(n)"
case1.err.test <- "rnorm(m)"

# produce testing set
m <- 100000 # sample size
x <- cbind(1, array(rnorm(m * d), c(m, d)))
ps <- exp(-0.5 - 0.5 * (x[, 2] + x[, 3])) / (1 + exp(-0.5 - 0.5 * (x[, 2] + x[, 3])))
a <- sapply(ps, rbinom, n = 1, size = 1)
eps <- eval(parse(text = case1.err.test))
y <- eval(parse(text = case1.base)) +   a * (x[, 2] + x[, 3]) + (2 * a + 1) * eps

# calculate the true optimal quantile 
case1.q.25 <- uniroot(opt.q.case1, c(-100, 100), tau = 0.25)$root
case1.q.50 <- uniroot(opt.q.case1, c(-100, 100), tau = 0.50)$root

# calculate the true optimal treatment regime for testing set
case1.opt.q.25 <- opt.reg.case1(case1.q.25)
case1.opt.q.50 <- opt.reg.case1(case1.q.50)

# parallel setup
cl <- makeCluster(getOption("cl.cores", clnum))
clusterSetRNGStream(cl, 222)
clusterEvalQ(cl, {
  library(MASS)
  library(glmnet)
  library(WeightSVM)
  library(quantoptr)
  library(quantreg)
})
clusterExport(cl, varlist = list("SCL.continuous", "SCL.continuous.hard", "QIQ.learning.continuous", "predict.QIQlearningcontinuous", "predict.glmnetformula"), envir = .GlobalEnv)

# parallel computing
system.time({
  res <- parLapply(cl, 1:times, simulation, n = n, d=d, err = case1.err, base = case1.base, base.test = case1.base.test, x.test = x, eps.test = eps, opt.q.25 = case1.opt.q.25, opt.q.50 = case1.opt.q.50, q.25 = case1.q.25, q.50 = case1.q.50)
})
stopCluster(cl)

# extract results
# value funtion and misclassification rates
case1.mianresult.SCL.linear <- array(0, dim = c(times, 4))
case1.mianresult.SCL.gaussian <- array(0, dim = c(times, 4))
case1.mianresult.QIQ <- array(0, dim = c(times, 4))
case1.mianresult.wang <- array(0, dim = c(times, 4))
case1.mianresult.SCL.linear.hard <- array(0, dim = c(times, 4))
case1.mianresult.SCL.gaussian.hard <- array(0, dim = c(times, 4))

#computation time
case1.time.SCL.linear <- array(0, dim = c(times, 2))
case1.time.SCL.gaussian <- array(0, dim = c(times, 2))
case1.time.QIQ <- array(0, dim = c(times, 2))
case1.time.wang <- array(0, dim = c(times, 2))

#bias of the estimated optimal value
case1.q.SCL.linear <- array(0, dim = c(times, 2))
case1.q.SCL.gaussian <- array(0, dim = c(times, 2))
case1.q.QIQ <- array(0, dim = c(times, 2))
case1.q.wang <- array(0, dim = c(times, 2))

#estimated value function
case1.hatv.SCL.linear <- array(0, dim = c(times, 2))
case1.hatv.SCL.gaussian <- array(0, dim = c(times, 2))
case1.hatv.QIQ <- array(0, dim = c(times, 2))
case1.hatv.wang <- array(0, dim = c(times, 2))

case1.hatv.true <- array(0, dim = c(times, 2))

if (d==2){
# quantile optimal treatment on a grid for plotting purpose
case1.xgrid.treat.SCL.linear.25 <- array(0, dim = c(times, 500*500))
case1.xgrid.treat.SCL.gaussian.25 <- array(0, dim = c(times, 500*500))
case1.xgrid.treat.QIQ.25 <- array(0, dim = c(times, 500*500))
case1.xgrid.treat.wang.25 <- array(0, dim = c(times, 500*500))

case1.xgrid.treat.SCL.linear.50 <- array(0, dim = c(times, 500*500))
case1.xgrid.treat.SCL.gaussian.50 <- array(0, dim = c(times, 500*500))
case1.xgrid.treat.QIQ.50 <- array(0, dim = c(times, 500*500))
case1.xgrid.treat.wang.50 <- array(0, dim = c(times, 500*500))
}

for (j in 1:times) {
  case1.mianresult.SCL.linear[j, ] <- res[[j]]$SCL.linear
  case1.mianresult.SCL.gaussian[j, ] <- res[[j]]$SCL.gaussian
  case1.mianresult.QIQ[j, ] <- res[[j]]$QIQ
  case1.mianresult.wang[j, ] <- res[[j]]$wang
  case1.mianresult.SCL.linear.hard[j, ] <- res[[j]]$SCL.linear.hard
  case1.mianresult.SCL.gaussian.hard[j, ] <- res[[j]]$SCL.gaussian.hard
   
  case1.time.SCL.linear[j, ] <- res[[j]]$time.SCL.linear
  case1.time.SCL.gaussian[j, ] <- res[[j]]$time.SCL.gaussian
  case1.time.QIQ[j, ] <- res[[j]]$time.QIQ
  case1.time.wang[j, ] <- res[[j]]$time.wang
 
  case1.q.SCL.linear[j, ] <- res[[j]]$q.SCL.linear
  case1.q.SCL.gaussian[j, ] <- res[[j]]$q.SCL.gaussian
  case1.q.QIQ[j, ] <- res[[j]]$q.QIQ
  case1.q.wang[j, ] <- res[[j]]$q.wang

  case1.hatv.SCL.linear[j, ] <- res[[j]]$hatv.SCL.linear
  case1.hatv.SCL.gaussian[j, ] <- res[[j]]$hatv.SCL.gaussian
  case1.hatv.QIQ[j, ] <- res[[j]]$hatv.QIQ
  case1.hatv.wang[j, ] <- res[[j]]$hatv.wang
  case1.hatv.true[j, ] <- res[[j]]$hatv.true
  
  if(d==2){
    case1.xgrid.treat.SCL.linear.25[j, ] <- res[[j]]$xgrid.treat.SCL.linear.25
    case1.xgrid.treat.SCL.gaussian.25[j, ] <- res[[j]]$xgrid.treat.SCL.gaussian.25
    case1.xgrid.treat.QIQ.25[j, ] <- res[[j]]$xgrid.treat.QIQ.25
    case1.xgrid.treat.wang.25[j, ] <- res[[j]]$xgrid.treat.wang.25

    case1.xgrid.treat.SCL.linear.50[j, ] <- res[[j]]$xgrid.treat.SCL.linear.50
    case1.xgrid.treat.SCL.gaussian.50[j, ] <- res[[j]]$xgrid.treat.SCL.gaussian.50
    case1.xgrid.treat.QIQ.50[j, ] <- res[[j]]$xgrid.treat.QIQ.50
    case1.xgrid.treat.wang.50[j, ] <- res[[j]]$xgrid.treat.wang.50
  }
}