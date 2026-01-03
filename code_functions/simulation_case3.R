# import necessary libraries
library(MASS)
library(parallel)
library(glmnet)
library(WeightSVM)
library(quantoptr)
library(mpath)

# source necessary functions
source("./code_functions/function_main.R")

### The main simulation function that generates data, fits different methods, and evaluates their performance.
### The input l denotes the index of the simulation run, the inputs-n, d, base-are used to generate data, and the inputs-base.test, x.test, opt.q.25, opt.q.50, q.25, q.50-are used to evaluate different methods. 
### The outputs include all metrics for competing methods, including the value, the estimated value, the bias of the estimated optimal value, misclassification rates, computation time, and quantile optimal treatment on a grid for plotting purposes
simulation <- function(l, n, d, base, base.test, x.test, opt.q.25, opt.q.50, q.25, q.50) {
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
  time.nuisance.ps.0 <- Sys.time()
  if(d<=5){
    f.a.x <- glmnet(x[, -1], a, family = binomial(link = "logit"), lambda=0)
    mean.a <- predict(f.a.x, x[, -1], type = "response")
  }else{
    f.a.x <- cv.glmnet(x[, -1], a, family = binomial(link = "logit"))
    mean.a <- predict(f.a.x, x[, -1], type = "response")
  }
  
  time.nuisance.ps.1 <- Sys.time()

  # estimation of the survival functions
  time.nuisance.sm.0 <- Sys.time()
  dataframe0 <- as.data.frame(list(y=y[a==0],x=x[a==0,-1]))
  dataframe1 <- as.data.frame(list(y=y[a==1],x=x[a==1,-1]))

  formula <-"y ~ "
  for(i in 1:d){
    formula = paste0(formula, paste("x.", i, sep = ""))
    if (i < d) {
      formula = paste0(formula, " + ")
    }
  }
  formula <- as.formula(formula)
  
  if(d<=5){
    f.y.ax0 <- glm.nb(formula, data = dataframe0)
    f.y.ax1 <- glm.nb(formula, data = dataframe1)
  }else{
    f.y.ax0 <- cv.glmregNB(formula, data = dataframe0, nfolds=5)
    f.y.ax1 <- cv.glmregNB(formula, data = dataframe1, nfolds=5)
    f.y.ax0$theta <- f.y.ax0$fit$theta[f.y.ax0$lambda.which]
    f.y.ax1$theta <- f.y.ax1$fit$theta[f.y.ax1$lambda.which]
  }

  time.nuisance.sm.1 <- Sys.time()

  ### Run competeting methods for estimating quantile OTRs
  # Our successive classification learning (SCL)
  kappa <- sd(y) / sqrt(n) / 6
  epsilon <- 1 / sqrt(n) / 2

  time_0 <- Sys.time()
  res.SCL.linear.25 <- SCL.count(data, 0.25, kappa, epsilon, mean.a, f.y.ax1, f.y.ax0, kernel="linear")
  time_1 <- Sys.time()
  res.SCL.gaussian.25 <- SCL.count(data, 0.25, kappa, epsilon, mean.a, f.y.ax1, f.y.ax0, kernel="gaussian")
  time_2 <- Sys.time()
  res.SCL.linear.50 <- SCL.count(data, 0.50, kappa, epsilon, mean.a, f.y.ax1, f.y.ax0, kernel="linear")
  time_3 <- Sys.time()
  res.SCL.gaussian.50 <- SCL.count(data, 0.50, kappa, epsilon, mean.a, f.y.ax1, f.y.ax0, kernel="gaussian")
  time_4 <- Sys.time()

  # QIQ learning
  res.QIQ.25 <- QIQ.learning.count(data, 0.25, f.y.ax1, f.y.ax0)
  time_5 <- Sys.time()
  res.QIQ.50 <- QIQ.learning.count(data, 0.50, f.y.ax1, f.y.ax0)
  time_6 <- Sys.time()

  # Wang's doubly robust method
  xformula <- ""
  for(j in 1:d){
    xformula <- paste(xformula,paste("x.", j, sep = ""),seq="")
    if (j < d) {
      xformula <- paste(xformula,"+",sep="")
    }
  }
  data.wang <- data.frame(y = y-runif(n), x = x[, -1], a = a)
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
  res.SCL.linear.25.hard <- SCL.count.hard(data, 0.25, kappa, epsilon, mean.a, f.y.ax1, f.y.ax0, kernel = "linear")
  res.SCL.gaussian.25.hard <- SCL.count.hard(data, 0.25, kappa, epsilon, mean.a, f.y.ax1, f.y.ax0, kernel = "gaussian")
  res.SCL.linear.50.hard <- SCL.count.hard(data, 0.50, kappa, epsilon, mean.a, f.y.ax1, f.y.ax0, kernel = "linear")
  res.SCL.gaussian.50.hard <- SCL.count.hard(data, 0.50, kappa, epsilon, mean.a, f.y.ax1, f.y.ax0, kernel = "gaussian")

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
  v.SCL.linear.25 <- smooth.quantile(rnbinom(100000, 2 * xtest.treat.SCL.linear.25 + 1, expit(eval(parse(text = base.test)) +  xtest.treat.SCL.linear.25 * (1 + x.test[, 2] + x.test[, 3]^2))), 0.25)
  v.SCL.linear.50 <- smooth.quantile(rnbinom(100000, 2 * xtest.treat.SCL.linear.50 + 1, expit(eval(parse(text = base.test)) +  xtest.treat.SCL.linear.50 * (1 + x.test[, 2] + x.test[, 3]^2))), 0.50)

  v.SCL.gaussian.25 <- smooth.quantile(rnbinom(100000, 2 * xtest.treat.SCL.gaussian.25 + 1, expit(eval(parse(text = base.test)) +  xtest.treat.SCL.gaussian.25 * (1 + x.test[, 2] + x.test[, 3]^2))), 0.25)
  v.SCL.gaussian.50 <- smooth.quantile(rnbinom(100000, 2 * xtest.treat.SCL.gaussian.50 + 1, expit(eval(parse(text = base.test)) +  xtest.treat.SCL.gaussian.50 * (1 + x.test[, 2] + x.test[, 3]^2))), 0.50)
  v.QIQ.25 <- smooth.quantile(rnbinom(100000, 2 * xtest.treat.QIQ.25 + 1, expit(eval(parse(text = base.test)) +  xtest.treat.QIQ.25 * (1 + x.test[, 2] + x.test[, 3]^2))), 0.25)
  v.QIQ.50 <- smooth.quantile(rnbinom(100000, 2 * xtest.treat.QIQ.50 + 1, expit(eval(parse(text = base.test)) +  xtest.treat.QIQ.50 * (1 + x.test[, 2] + x.test[, 3]^2))), 0.50)

  v.wang.25 <- smooth.quantile(rnbinom(100000, 2 * xtest.treat.wang.25 + 1, expit(eval(parse(text = base.test)) +  xtest.treat.wang.25 * (1 + x.test[, 2] + x.test[, 3]^2))), 0.25)
  v.wang.50 <- smooth.quantile(rnbinom(100000, 2 * xtest.treat.wang.50 + 1, expit(eval(parse(text = base.test)) +  xtest.treat.wang.50 * (1 + x.test[, 2] + x.test[, 3]^2))), 0.50)

  v.SCL.linear.25.hard <- smooth.quantile(rnbinom(100000, 2 * xtest.treat.SCL.linear.25.hard + 1, expit(eval(parse(text = base.test)) +  xtest.treat.SCL.linear.25.hard * (1 + x.test[, 2] + x.test[, 3]^2))), 0.25)
  v.SCL.linear.50.hard <- smooth.quantile(rnbinom(100000, 2 * xtest.treat.SCL.linear.50.hard + 1, expit(eval(parse(text = base.test)) +  xtest.treat.SCL.linear.50.hard * (1 + x.test[, 2] + x.test[, 3]^2))), 0.50)
  v.SCL.gaussian.25.hard <- smooth.quantile(rnbinom(100000, 2 * xtest.treat.SCL.gaussian.25.hard + 1, expit(eval(parse(text = base.test)) +  xtest.treat.SCL.gaussian.25.hard * (1 + x.test[, 2] + x.test[, 3]^2))), 0.25)
  v.SCL.gaussian.50.hard <- smooth.quantile(rnbinom(100000, 2 * xtest.treat.SCL.gaussian.50.hard + 1, expit(eval(parse(text = base.test)) +  xtest.treat.SCL.gaussian.50.hard * (1 + x.test[, 2] + x.test[, 3]^2))), 0.50)

  ## estimated value function
  hatv.SCL.linear.25 <-  hat.smooth.quantile(y, a, mean.a, xtreat.SCL.linear.25, 0.25)
  hatv.SCL.linear.50 <-  hat.smooth.quantile(y, a, mean.a, xtreat.SCL.linear.50, 0.50)

  hatv.SCL.gaussian.25 <-  hat.smooth.quantile(y, a, mean.a, xtreat.SCL.gaussian.25, 0.25)
  hatv.SCL.gaussian.50 <-  hat.smooth.quantile(y, a, mean.a, xtreat.SCL.gaussian.50, 0.50)

  hatv.QIQ.25 <-  hat.smooth.quantile(y, a, mean.a, xtreat.QIQ.25, 0.25)
  hatv.QIQ.50 <-  hat.smooth.quantile(y, a, mean.a, xtreat.QIQ.50, 0.50)

  hatv.wang.25 <-  hat.smooth.quantile(y, a, mean.a, xtreat.wang.25, 0.25)
  hatv.wang.50 <-  hat.smooth.quantile(y, a, mean.a, xtreat.wang.50, 0.50)
  
  # calculate the true quantile treatment for x
  qfloor <- floor(q.25)
  qceiling <- ceiling(q.25)

  if (qfloor == qceiling) {
    lambda <- 0
  } else {
    lambda <- (q.25 - qceiling) / (qfloor - qceiling)
  }

  S0 <- pnbinom(qceiling, prob = expit(eval(parse(text = base))), lower.tail = FALSE, size = 1) + lambda * dnbinom(qceiling, prob = expit(eval(parse(text = base))), size = 1)

  S1 <- pnbinom(qceiling, prob = expit(eval(parse(text = base)) + 1 + x[, 2] + x[, 3]^2), lower.tail = FALSE, size = 3) + lambda * dnbinom(qceiling, prob = expit(eval(parse(text = base)) + 1 + x[, 2] + x[, 3]^2), size = 3)

  true.25 <- as.numeric((S1-S0) > 0)

  qfloor <- floor(q.50)
  qceiling <- ceiling(q.50)

  if (qfloor == qceiling) {
    lambda <- 0
  } else {
    lambda <- (q.50 - qceiling) / (qfloor - qceiling)
  }

  S0 <- pnbinom(qceiling, prob = expit(eval(parse(text = base))), lower.tail = FALSE, size = 1) + lambda * dnbinom(qceiling, prob = expit(eval(parse(text = base))), size = 1)

  S1 <- pnbinom(qceiling, prob = expit(eval(parse(text = base)) + 1 + x[, 2] + x[, 3]^2), lower.tail = FALSE, size = 3) + lambda * dnbinom(qceiling, prob = expit(eval(parse(text = base)) + 1+ x[, 2] + x[, 3]^2), size = 3)

  true.50 <- as.numeric((S1-S0) > 0)

  # estimated value function based on the true quantile treatment
  hatv.true.25 <- hat.smooth.quantile(y, a, mean.a, true.25, 0.25)
  hatv.true.50 <- hat.smooth.quantile(y, a, mean.a, true.50, 0.50)

  # Calculate quantile optimal treatment on a grid for plotting purpose, which is only feasible when d=2.
  if(d==2){
    df <- expand.grid(x1 = seq(-2, 2, length.out = 500),
                    x2 = seq(-2, 2, length.out = 500))
    x.grid <- cbind(1,df,deparse.level = 0)
    names(x.grid)[2] <- "x.1"
    names(x.grid)[3] <- "x.2"

    xgrid.treat.SCL.linear.25 <- as.numeric(predict(res.SCL.linear.25$opt.reg, x.grid[,-1])) - 1
    xgrid.treat.SCL.gaussian.25 <- as.numeric(predict(res.SCL.gaussian.25$opt.reg, x.grid[,-1])) - 1
    xgrid.treat.SCL.linear.50 <- as.numeric(predict(res.SCL.linear.50$opt.reg, x.grid[,-1])) - 1
    xgrid.treat.SCL.gaussian.50 <- as.numeric(predict(res.SCL.gaussian.50$opt.reg, x.grid[,-1])) - 1
    xgrid.treat.QIQ.25 <- as.numeric(predict(res.QIQ.25$opt.reg, unname(x.grid[,-1])))
    xgrid.treat.QIQ.50 <- as.numeric(predict(res.QIQ.50$opt.reg, unname(x.grid[,-1])))

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
      xgrid.treat.wang.50= xgrid.treat.wang.50
    ))} else {
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
    ))}
}

# Specifiy the data generating process. If you want to test other data generating processes, please change the base, n accordingly.
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

# calculate the true optimal quantile 
case3.q.25 <- uniroot(opt.q.case3, c(-100, 100), tau = 0.25)$root
case3.q.50 <- uniroot(opt.q.case3, c(-100, 100), tau = 0.50)$root

# calculate the true optimal treatment regime for testing set
case3.opt.q.25 <- opt.reg.case3(case3.q.25)
case3.opt.q.50 <- opt.reg.case3(case3.q.50)

# parallel setup
cl <- makeCluster(getOption("cl.cores", clnum))
clusterSetRNGStream(cl, 222)
clusterEvalQ(cl, {
  library(MASS)
  library(glmnet)
  library(WeightSVM)
  library(quantoptr)
  library(mpath)
})
clusterExport(cl, varlist = list("SCL.count", "SCL.count.hard", "QIQ.learning.count", "predict.QIQlearningcount", "expit", "smooth.quantile", "hat.smooth.quantile"), envir = .GlobalEnv)

# parallel computing
system.time({
  res <- parLapply(cl, 1:times, simulation, n = n, d=d, base = case3.base, base.test = case3.base.test, x.test = x, opt.q.25 = case3.opt.q.25, opt.q.50 = case3.opt.q.50, q.25 = case3.q.25,q.50 = case3.q.50)
})
stopCluster(cl)

# extract results
# value function and misclassification rates
case3.mianresult.SCL.linear <- array(0, dim = c(times, 4))
case3.mianresult.SCL.gaussian <- array(0, dim = c(times, 4))
case3.mianresult.QIQ <- array(0, dim = c(times, 4))
case3.mianresult.wang <- array(0, dim = c(times, 4))
case3.mianresult.SCL.linear.hard <- array(0, dim = c(times, 4))
case3.mianresult.SCL.gaussian.hard <- array(0, dim = c(times, 4))

#computation time
case3.time.SCL.linear <- array(0, dim = c(times, 2))
case3.time.SCL.gaussian <- array(0, dim = c(times, 2))
case3.time.QIQ <- array(0, dim = c(times, 2))
case3.time.wang <- array(0, dim = c(times, 2))

#bias of the estimated optimal value
case3.q.SCL.linear <- array(0, dim = c(times, 2))
case3.q.SCL.gaussian <- array(0, dim = c(times, 2))
case3.q.QIQ <- array(0, dim = c(times, 2))
case3.q.wang <- array(0, dim = c(times, 2))

#estimated value function
case3.hatv.SCL.linear <- array(0, dim = c(times, 2))
case3.hatv.SCL.gaussian <- array(0, dim = c(times, 2))
case3.hatv.QIQ <- array(0, dim = c(times, 2))
case3.hatv.wang <- array(0, dim = c(times, 2))

case3.hatv.true <- array(0, dim = c(times, 2))

# quantile optimal treatment on a grid for plotting purpose, which is only feasible when d=2.
if(d==2){
  case3.xgrid.treat.SCL.linear.25 <- array(0, dim = c(times, 500*500))
  case3.xgrid.treat.SCL.gaussian.25 <- array(0, dim = c(times, 500*500))
  case3.xgrid.treat.QIQ.25 <- array(0, dim = c(times, 500*500))
  case3.xgrid.treat.wang.25 <- array(0, dim = c(times, 500*500))

  case3.xgrid.treat.SCL.linear.50 <- array(0, dim = c(times, 500*500))
  case3.xgrid.treat.SCL.gaussian.50 <- array(0, dim = c(times, 500*500))
  case3.xgrid.treat.QIQ.50 <- array(0, dim = c(times, 500*500))
  case3.xgrid.treat.wang.50 <- array(0, dim = c(times, 500*500))
}
for (j in 1:times) {
  case3.mianresult.SCL.linear[j, ] <- res[[j]]$SCL.linear
  case3.mianresult.SCL.gaussian[j, ] <- res[[j]]$SCL.gaussian
  case3.mianresult.QIQ[j, ] <- res[[j]]$QIQ
  case3.mianresult.wang[j, ] <- res[[j]]$wang
  case3.mianresult.SCL.linear.hard[j, ] <- res[[j]]$SCL.linear.hard
  case3.mianresult.SCL.gaussian.hard[j, ] <- res[[j]]$SCL.gaussian.hard
  
  case3.time.SCL.linear[j, ] <- res[[j]]$time.SCL.linear
  case3.time.SCL.gaussian[j, ] <- res[[j]]$time.SCL.gaussian
  case3.time.QIQ[j, ] <- res[[j]]$time.QIQ
  case3.time.wang[j, ] <- res[[j]]$time.wang
 
  case3.q.SCL.linear[j, ] <- res[[j]]$q.SCL.linear
  case3.q.SCL.gaussian[j, ] <- res[[j]]$q.SCL.gaussian
  case3.q.QIQ[j, ] <- res[[j]]$q.QIQ
  case3.q.wang[j, ] <- res[[j]]$q.wang

  case3.hatv.SCL.linear[j, ] <- res[[j]]$hatv.SCL.linear
  case3.hatv.SCL.gaussian[j, ] <- res[[j]]$hatv.SCL.gaussian
  case3.hatv.QIQ[j, ] <- res[[j]]$hatv.QIQ
  case3.hatv.wang[j, ] <- res[[j]]$hatv.wang
  case3.hatv.true[j, ] <- res[[j]]$hatv.true
  
  if(d==2){
    case3.xgrid.treat.SCL.linear.25[j, ] <- res[[j]]$xgrid.treat.SCL.linear.25
    case3.xgrid.treat.SCL.gaussian.25[j, ] <- res[[j]]$xgrid.treat.SCL.gaussian.25
    case3.xgrid.treat.QIQ.25[j, ] <- res[[j]]$xgrid.treat.QIQ.25
    case3.xgrid.treat.wang.25[j, ] <- res[[j]]$xgrid.treat.wang.25

    case3.xgrid.treat.SCL.linear.50[j, ] <- res[[j]]$xgrid.treat.SCL.linear.50
    case3.xgrid.treat.SCL.gaussian.50[j, ] <- res[[j]]$xgrid.treat.SCL.gaussian.50
    case3.xgrid.treat.QIQ.50[j, ] <- res[[j]]$xgrid.treat.QIQ.50
    case3.xgrid.treat.wang.50[j, ] <- res[[j]]$xgrid.treat.wang.50
  }
}