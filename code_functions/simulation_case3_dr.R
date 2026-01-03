library(MASS)
library(parallel)
library(glmnet)
library(WeightSVM)
library(quantoptr)

# source necessary functions
source("./code_functions/function_main.R")

### The main simulation function that generates data, fits different methods, and evaluates their performance.
### The input l denotes the index of the simulation run, the inputs-n, d, base-are used to generate data, and the inputs-base.test, x.test, eps.test, opt.q.25, opt.q.50, q.25, q.50-are used to evaluate different methods. 
### The outputs include the value and the misclassification rate for different methods.
simulation <- function(l, n, d, base, base.test, x.test, eps.test, opt.q.25, opt.q.50) {
  set.seed(l)
  # generate data
  n <- n
  d <- d

  x <- cbind(1, array(rnorm(n * d), c(n, d)))

  ps <- exp(-0.5 - 0.5 * (x[,2]+x[,3])) / (1 + exp(-0.5 - 0.5 * (x[,2]+x[,3])))
  a <- sapply(ps, rbinom, n = 1, size = 1)
  
  mean.y <- expit(eval(parse(text = base)) +  a * (1 + x[, 2] + x[, 3]^2))
  y <- rnbinom(n, size = 2 * a + 1, prob = mean.y)

  data <- list(y = y, x = x, a = a)

  # estimation of the propensity score model
  f.a.x.correct <- glmnet(x[, -1], a, family = binomial(link = "logit"), lambda=0)
  mean.a.correct <- predict(f.a.x.correct, x[, -1], type = "response")

  # estimation of the incorrect propensity score model1
  f.a.x.incorrect1 <- glm(a ~ 1, family = binomial(link = "logit"))
  mean.a.incorrect1 <- predict(f.a.x.incorrect1, type = "response")

  # estimation of the incorrect propensity score model2
  constructx <- cbind(exp(x[,2]), x[, 3])
  f.a.x.incorrect2 <- glm(a ~ constructx, family = binomial(link = "logit"))
  mean.a.incorrect2 <- predict(f.a.x.incorrect2, type = "response")

  # estimation of the correct conditional survival model
  dataframe0 <- as.data.frame(list(y=y[a==0],x=x[a==0,-1]))
  dataframe1 <- as.data.frame(list(y=y[a==1],x=x[a==1,-1]))

  formula <- as.formula("y ~ x.1 + x.2 +I(x.2^2)")
  f.y.ax0.correct <- glm.nb(formula, data = dataframe0)
  f.y.ax1.correct <- glm.nb(formula, data = dataframe1)
  
  # estimation of the incorrect conditional survival model
  formula <- as.formula("y ~ x.1 + x.2")
  f.y.ax0.incorrect <- glm.nb(formula, data = dataframe0)
  f.y.ax1.incorrect <- glm.nb(formula, data = dataframe1)
  
  # successive classification learning
  kappa <- sd(y) / sqrt(n) / 6
  epsilon <- 1 / sqrt(n) / 2

  # correct ps correct sm
  res.SCL.linear.25.correctcorrect <- SCL.count(data, 0.25, kappa, epsilon, mean.a.correct, f.y.ax1.correct, f.y.ax0.correct, kernel="linear")
  res.SCL.gaussian.25.correctcorrect <- SCL.count(data, 0.25, kappa, epsilon, mean.a.correct, f.y.ax1.correct, f.y.ax0.correct, kernel="gaussian")
  res.SCL.linear.50.correctcorrect <- SCL.count(data, 0.50, kappa, epsilon, mean.a.correct, f.y.ax1.correct, f.y.ax0.correct, kernel="linear")
  res.SCL.gaussian.50.correctcorrect <- SCL.count(data, 0.50, kappa, epsilon, mean.a.correct, f.y.ax1.correct, f.y.ax0.correct, kernel="gaussian")

  # incorrect ps1 correct sm
  res.SCL.linear.25.incorrect1correct <- SCL.count(data, 0.25, kappa, epsilon, mean.a.incorrect1, f.y.ax1.correct, f.y.ax0.correct, kernel="linear")
  res.SCL.gaussian.25.incorrect1correct <- SCL.count(data, 0.25, kappa, epsilon, mean.a.incorrect1, f.y.ax1.correct, f.y.ax0.correct, kernel="gaussian")
  res.SCL.linear.50.incorrect1correct <- SCL.count(data, 0.50, kappa, epsilon, mean.a.incorrect1, f.y.ax1.correct, f.y.ax0.correct, kernel="linear")
  res.SCL.gaussian.50.incorrect1correct <- SCL.count(data, 0.50, kappa, epsilon, mean.a.incorrect1, f.y.ax1.correct, f.y.ax0.correct, kernel="gaussian")

  # incorrect ps2 correct sm
  res.SCL.linear.25.incorrect2correct <- SCL.count(data, 0.25, kappa, epsilon, mean.a.incorrect2, f.y.ax1.correct, f.y.ax0.correct, kernel="linear")
  res.SCL.gaussian.25.incorrect2correct <- SCL.count(data, 0.25, kappa, epsilon, mean.a.incorrect2, f.y.ax1.correct, f.y.ax0.correct, kernel="gaussian")
  res.SCL.linear.50.incorrect2correct <- SCL.count(data, 0.50, kappa, epsilon, mean.a.incorrect2, f.y.ax1.correct, f.y.ax0.correct, kernel="linear")
  res.SCL.gaussian.50.incorrect2correct <- SCL.count(data, 0.50, kappa, epsilon, mean.a.incorrect2, f.y.ax1.correct, f.y.ax0.correct, kernel="gaussian")

  # correct ps incorrect sm
  res.SCL.linear.25.correctincorrect <- SCL.count(data, 0.25, kappa, epsilon, mean.a.correct, f.y.ax1.incorrect, f.y.ax0.incorrect, kernel="linear")
  res.SCL.gaussian.25.correctincorrect <- SCL.count(data, 0.25, kappa, epsilon, mean.a.correct, f.y.ax1.incorrect, f.y.ax0.incorrect, kernel="gaussian")
  res.SCL.linear.50.correctincorrect <- SCL.count(data, 0.50, kappa, epsilon, mean.a.correct, f.y.ax1.incorrect, f.y.ax0.incorrect, kernel="linear")
  res.SCL.gaussian.50.correctincorrect <- SCL.count(data, 0.50, kappa, epsilon, mean.a.correct, f.y.ax1.incorrect, f.y.ax0.incorrect, kernel="gaussian")

  # incorrect ps1 incorrect sm
  res.SCL.linear.25.incorrect1incorrect <- SCL.count(data, 0.25, kappa, epsilon, mean.a.incorrect1, f.y.ax1.incorrect, f.y.ax0.incorrect, kernel="linear")
  res.SCL.gaussian.25.incorrect1incorrect <- SCL.count(data, 0.25, kappa, epsilon, mean.a.incorrect1, f.y.ax1.incorrect, f.y.ax0.incorrect, kernel="gaussian")
  res.SCL.linear.50.incorrect1incorrect <- SCL.count(data, 0.50, kappa, epsilon, mean.a.incorrect1, f.y.ax1.incorrect, f.y.ax0.incorrect, kernel="linear")
  res.SCL.gaussian.50.incorrect1incorrect <- SCL.count(data, 0.50, kappa, epsilon, mean.a.incorrect1, f.y.ax1.incorrect, f.y.ax0.incorrect, kernel="gaussian")

  # incorrect ps2 incorrect sm
  res.SCL.linear.25.incorrect2incorrect <- SCL.count(data, 0.25, kappa, epsilon, mean.a.incorrect2, f.y.ax1.incorrect, f.y.ax0.incorrect, kernel="linear")
  res.SCL.gaussian.25.incorrect2incorrect <- SCL.count(data, 0.25, kappa, epsilon, mean.a.incorrect2, f.y.ax1.incorrect, f.y.ax0.incorrect, kernel="gaussian")
  res.SCL.linear.50.incorrect2incorrect <- SCL.count(data, 0.50, kappa, epsilon, mean.a.incorrect2, f.y.ax1.incorrect, f.y.ax0.incorrect, kernel="linear")
  res.SCL.gaussian.50.incorrect2incorrect <- SCL.count(data, 0.50, kappa, epsilon, mean.a.incorrect2, f.y.ax1.incorrect, f.y.ax0.incorrect, kernel="gaussian")

  #extract result for prediction
  treat.SCL.linear.25.correctcorrect <- as.numeric(predict(res.SCL.linear.25.correctcorrect$opt.reg, x.test[, -1])) - 1
  treat.SCL.gaussian.25.correctcorrect <- as.numeric(predict(res.SCL.gaussian.25.correctcorrect$opt.reg, x.test[, -1])) - 1
  treat.SCL.linear.50.correctcorrect <- as.numeric(predict(res.SCL.linear.50.correctcorrect$opt.reg, x.test[, -1])) - 1
  treat.SCL.gaussian.50.correctcorrect <- as.numeric(predict(res.SCL.gaussian.50.correctcorrect$opt.reg, x.test[, -1])) - 1

  treat.SCL.linear.25.incorrect1correct <- as.numeric(predict(res.SCL.linear.25.incorrect1correct$opt.reg, x.test[, -1])) - 1
  treat.SCL.gaussian.25.incorrect1correct <- as.numeric(predict(res.SCL.gaussian.25.incorrect1correct$opt.reg, x.test[, -1])) - 1
  treat.SCL.linear.50.incorrect1correct <- as.numeric(predict(res.SCL.linear.50.incorrect1correct$opt.reg, x.test[, -1])) - 1
  treat.SCL.gaussian.50.incorrect1correct <- as.numeric(predict(res.SCL.gaussian.50.incorrect1correct$opt.reg, x.test[, -1])) - 1

  treat.SCL.linear.25.incorrect2correct <- as.numeric(predict(res.SCL.linear.25.incorrect2correct$opt.reg, x.test[, -1])) - 1
  treat.SCL.gaussian.25.incorrect2correct <- as.numeric(predict(res.SCL.gaussian.25.incorrect2correct$opt.reg, x.test[, -1])) - 1
  treat.SCL.linear.50.incorrect2correct <- as.numeric(predict(res.SCL.linear.50.incorrect2correct$opt.reg, x.test[, -1])) - 1
  treat.SCL.gaussian.50.incorrect2correct <- as.numeric(predict(res.SCL.gaussian.50.incorrect2correct$opt.reg, x.test[, -1])) - 1

  treat.SCL.linear.25.correctincorrect <- as.numeric(predict(res.SCL.linear.25.correctincorrect$opt.reg, x.test[, -1])) - 1
  treat.SCL.gaussian.25.correctincorrect <- as.numeric(predict(res.SCL.gaussian.25.correctincorrect$opt.reg, x.test[, -1])) - 1
  treat.SCL.linear.50.correctincorrect <- as.numeric(predict(res.SCL.linear.50.correctincorrect$opt.reg, x.test[, -1])) - 1
  treat.SCL.gaussian.50.correctincorrect <- as.numeric(predict(res.SCL.gaussian.50.correctincorrect$opt.reg, x.test[, -1])) - 1

  treat.SCL.linear.25.incorrect1incorrect <- as.numeric(predict(res.SCL.linear.25.incorrect1incorrect$opt.reg, x.test[, -1])) - 1
  treat.SCL.gaussian.25.incorrect1incorrect <- as.numeric(predict(res.SCL.gaussian.25.incorrect1incorrect$opt.reg, x.test[, -1])) - 1
  treat.SCL.linear.50.incorrect1incorrect <- as.numeric(predict(res.SCL.linear.50.incorrect1incorrect$opt.reg, x.test[, -1])) - 1
  treat.SCL.gaussian.50.incorrect1incorrect <- as.numeric(predict(res.SCL.gaussian.50.incorrect1incorrect$opt.reg, x.test[, -1])) - 1

  treat.SCL.linear.25.incorrect2incorrect <- as.numeric(predict(res.SCL.linear.25.incorrect2incorrect$opt.reg, x.test[, -1])) - 1
  treat.SCL.gaussian.25.incorrect2incorrect <- as.numeric(predict(res.SCL.gaussian.25.incorrect2incorrect$opt.reg, x.test[, -1])) - 1
  treat.SCL.linear.50.incorrect2incorrect <- as.numeric(predict(res.SCL.linear.50.incorrect2incorrect$opt.reg, x.test[, -1])) - 1
  treat.SCL.gaussian.50.incorrect2incorrect <- as.numeric(predict(res.SCL.gaussian.50.incorrect2incorrect$opt.reg, x.test[, -1])) - 1

  # misclassification rate
  mr.SCL.linear.25.correctcorrect <- 1 - sum(treat.SCL.linear.25.correctcorrect == opt.q.25) / dim(x.test)[1]
  mr.SCL.gaussian.25.correctcorrect <- 1 - sum(treat.SCL.gaussian.25.correctcorrect == opt.q.25) / dim(x.test)[1]
  mr.SCL.linear.50.correctcorrect <- 1 - sum(treat.SCL.linear.50.correctcorrect == opt.q.50) / dim(x.test)[1]
  mr.SCL.gaussian.50.correctcorrect <- 1 - sum(treat.SCL.gaussian.50.correctcorrect == opt.q.50) / dim(x.test)[1]

  mr.SCL.linear.25.incorrect1correct <- 1 - sum(treat.SCL.linear.25.incorrect1correct == opt.q.25) / dim(x.test)[1]
  mr.SCL.gaussian.25.incorrect1correct <- 1 - sum(treat.SCL.gaussian.25.incorrect1correct == opt.q.25) / dim(x.test)[1]
  mr.SCL.linear.50.incorrect1correct <- 1 - sum(treat.SCL.linear.50.incorrect1correct == opt.q.50) / dim(x.test)[1]
  mr.SCL.gaussian.50.incorrect1correct <- 1 - sum(treat.SCL.gaussian.50.incorrect1correct == opt.q.50) / dim(x.test)[1]

  mr.SCL.linear.25.incorrect2correct <- 1 - sum(treat.SCL.linear.25.incorrect2correct == opt.q.25) / dim(x.test)[1]
  mr.SCL.gaussian.25.incorrect2correct <- 1 - sum(treat.SCL.gaussian.25.incorrect2correct == opt.q.25) / dim(x.test)[1]
  mr.SCL.linear.50.incorrect2correct <- 1 - sum(treat.SCL.linear.50.incorrect2correct == opt.q.50) / dim(x.test)[1]
  mr.SCL.gaussian.50.incorrect2correct <- 1 - sum(treat.SCL.gaussian.50.incorrect2correct == opt.q.50) / dim(x.test)[1]

  mr.SCL.linear.25.correctincorrect <- 1 - sum(treat.SCL.linear.25.correctincorrect == opt.q.25) / dim(x.test)[1]
  mr.SCL.gaussian.25.correctincorrect <- 1 - sum(treat.SCL.gaussian.25.correctincorrect == opt.q.25) / dim(x.test)[1]
  mr.SCL.linear.50.correctincorrect <- 1 - sum(treat.SCL.linear.50.correctincorrect == opt.q.50) / dim(x.test)[1]
  mr.SCL.gaussian.50.correctincorrect <- 1 - sum(treat.SCL.gaussian.50.correctincorrect == opt.q.50) / dim(x.test)[1]

  mr.SCL.linear.25.incorrect1incorrect <- 1 - sum(treat.SCL.linear.25.incorrect1incorrect == opt.q.25) / dim(x.test)[1]
  mr.SCL.gaussian.25.incorrect1incorrect <- 1 - sum(treat.SCL.gaussian.25.incorrect1incorrect == opt.q.25) / dim(x.test)[1]
  mr.SCL.linear.50.incorrect1incorrect <- 1 - sum(treat.SCL.linear.50.incorrect1incorrect == opt.q.50) / dim(x.test)[1]
  mr.SCL.gaussian.50.incorrect1incorrect <- 1 - sum(treat.SCL.gaussian.50.incorrect1incorrect == opt.q.50) / dim(x.test)[1]

  mr.SCL.linear.25.incorrect2incorrect <- 1 - sum(treat.SCL.linear.25.incorrect2incorrect == opt.q.25) / dim(x.test)[1]
  mr.SCL.gaussian.25.incorrect2incorrect <- 1 - sum(treat.SCL.gaussian.25.incorrect2incorrect == opt.q.25) / dim(x.test)[1]
  mr.SCL.linear.50.incorrect2incorrect <- 1 - sum(treat.SCL.linear.50.incorrect2incorrect == opt.q.50) / dim(x.test)[1]
  mr.SCL.gaussian.50.incorrect2incorrect <- 1 - sum(treat.SCL.gaussian.50.incorrect2incorrect == opt.q.50) / dim(x.test)[1]

  # value function
  set.seed(l)
  v.SCL.linear.25.correctcorrect <- smooth.quantile(rnbinom(100000, 2 * treat.SCL.linear.25.correctcorrect + 1, expit(eval(parse(text = base.test)) +  treat.SCL.linear.25.correctcorrect * (1 + x.test[, 2] + x.test[, 3]^2))), 0.25)
  v.SCL.linear.50.correctcorrect <- smooth.quantile(rnbinom(100000, 2 * treat.SCL.linear.50.correctcorrect + 1, expit(eval(parse(text = base.test)) +  treat.SCL.linear.50.correctcorrect * (1 + x.test[, 2] + x.test[, 3]^2))), 0.50)
  v.SCL.gaussian.25.correctcorrect <- smooth.quantile(rnbinom(100000, 2 * treat.SCL.gaussian.25.correctcorrect + 1, expit(eval(parse(text = base.test)) +  treat.SCL.gaussian.25.correctcorrect * (1 + x.test[, 2] + x.test[, 3]^2))), 0.25)
  v.SCL.gaussian.50.correctcorrect <- smooth.quantile(rnbinom(100000, 2 * treat.SCL.gaussian.50.correctcorrect + 1, expit(eval(parse(text = base.test)) +  treat.SCL.gaussian.50.correctcorrect * (1 + x.test[, 2] + x.test[, 3]^2))), 0.50)

  v.SCL.linear.25.incorrect1correct <- smooth.quantile(rnbinom(100000, 2 * treat.SCL.linear.25.incorrect1correct + 1, expit(eval(parse(text = base.test)) +  treat.SCL.linear.25.incorrect1correct * (1 + x.test[, 2] + x.test[, 3]^2))), 0.25)
  v.SCL.gaussian.25.incorrect1correct <- smooth.quantile(rnbinom(100000, 2 * treat.SCL.gaussian.25.incorrect1correct + 1, expit(eval(parse(text = base.test)) +  treat.SCL.gaussian.25.incorrect1correct * (1 + x.test[, 2] + x.test[, 3]^2))), 0.25)
  v.SCL.linear.50.incorrect1correct <- smooth.quantile(rnbinom(100000, 2 * treat.SCL.linear.50.incorrect1correct + 1, expit(eval(parse(text = base.test)) +  treat.SCL.linear.50.incorrect1correct * (1 + x.test[, 2] + x.test[, 3]^2))), 0.50)
  v.SCL.gaussian.50.incorrect1correct <- smooth.quantile(rnbinom(100000, 2 * treat.SCL.gaussian.50.incorrect1correct + 1, expit(eval(parse(text = base.test)) +  treat.SCL.gaussian.50.incorrect1correct * (1 + x.test[, 2] + x.test[, 3]^2))), 0.50)

  v.SCL.linear.25.incorrect2correct <- smooth.quantile(rnbinom(100000, 2 * treat.SCL.linear.25.incorrect2correct + 1, expit(eval(parse(text = base.test)) +  treat.SCL.linear.25.incorrect2correct * (1 + x.test[, 2] + x.test[, 3]^2))), 0.25)
  v.SCL.gaussian.25.incorrect2correct <- smooth.quantile(rnbinom(100000, 2 * treat.SCL.gaussian.25.incorrect2correct + 1, expit(eval(parse(text = base.test)) +  treat.SCL.gaussian.25.incorrect2correct * (1 + x.test[, 2] + x.test[, 3]^2))), 0.25)
  v.SCL.linear.50.incorrect2correct <- smooth.quantile(rnbinom(100000, 2 * treat.SCL.linear.50.incorrect2correct + 1, expit(eval(parse(text = base.test)) +  treat.SCL.linear.50.incorrect2correct * (1 + x.test[, 2] + x.test[, 3]^2))), 0.50)
  v.SCL.gaussian.50.incorrect2correct <- smooth.quantile(rnbinom(100000, 2 * treat.SCL.gaussian.50.incorrect2correct + 1, expit(eval(parse(text = base.test)) +  treat.SCL.gaussian.50.incorrect2correct * (1 + x.test[, 2] + x.test[, 3]^2))), 0.50)

  v.SCL.linear.25.correctincorrect <- smooth.quantile(rnbinom(100000, 2 * treat.SCL.linear.25.correctincorrect + 1, expit(eval(parse(text = base.test)) +  treat.SCL.linear.25.correctincorrect * (1 + x.test[, 2] + x.test[, 3]^2))), 0.25)
  v.SCL.gaussian.25.correctincorrect <- smooth.quantile(rnbinom(100000, 2 * treat.SCL.gaussian.25.correctincorrect + 1, expit(eval(parse(text = base.test)) +  treat.SCL.gaussian.25.correctincorrect * (1 + x.test[, 2] + x.test[, 3]^2))), 0.25)
  v.SCL.linear.50.correctincorrect <- smooth.quantile(rnbinom(100000, 2 * treat.SCL.linear.50.correctincorrect + 1, expit(eval(parse(text = base.test)) +  treat.SCL.linear.50.correctincorrect * (1 + x.test[, 2] + x.test[, 3]^2))), 0.50)
  v.SCL.gaussian.50.correctincorrect <- smooth.quantile(rnbinom(100000, 2 * treat.SCL.gaussian.50.correctincorrect + 1, expit(eval(parse(text = base.test)) +  treat.SCL.gaussian.50.correctincorrect * (1 + x.test[, 2] + x.test[, 3]^2))), 0.50)

  v.SCL.linear.25.incorrect1incorrect <- smooth.quantile(rnbinom(100000, 2 * treat.SCL.linear.25.incorrect1incorrect + 1, expit(eval(parse(text = base.test)) +  treat.SCL.linear.25.incorrect1incorrect * (1 + x.test[, 2] + x.test[, 3]^2))), 0.25)
  v.SCL.gaussian.25.incorrect1incorrect <- smooth.quantile(rnbinom(100000, 2 * treat.SCL.gaussian.25.incorrect1incorrect + 1, expit(eval(parse(text = base.test)) +  treat.SCL.gaussian.25.incorrect1incorrect * (1 + x.test[, 2] + x.test[, 3]^2))), 0.25)
  v.SCL.linear.50.incorrect1incorrect <- smooth.quantile(rnbinom(100000, 2 * treat.SCL.linear.50.incorrect1incorrect + 1, expit(eval(parse(text = base.test)) +  treat.SCL.linear.50.incorrect1incorrect * (1 + x.test[, 2] + x.test[, 3]^2))), 0.50)
  v.SCL.gaussian.50.incorrect1incorrect <- smooth.quantile(rnbinom(100000, 2 * treat.SCL.gaussian.50.incorrect1incorrect + 1, expit(eval(parse(text = base.test)) +  treat.SCL.gaussian.50.incorrect1incorrect * (1 + x.test[, 2] + x.test[, 3]^2))), 0.50)

  v.SCL.linear.25.incorrect2incorrect <- smooth.quantile(rnbinom(100000, 2 * treat.SCL.linear.25.incorrect2incorrect + 1, expit(eval(parse(text = base.test)) +  treat.SCL.linear.25.incorrect2incorrect * (1 + x.test[, 2] + x.test[, 3]^2))), 0.25)
  v.SCL.gaussian.25.incorrect2incorrect <- smooth.quantile(rnbinom(100000, 2 * treat.SCL.gaussian.25.incorrect2incorrect + 1, expit(eval(parse(text = base.test)) +  treat.SCL.gaussian.25.incorrect2incorrect * (1 + x.test[, 2] + x.test[, 3]^2))), 0.25)
  v.SCL.linear.50.incorrect2incorrect <- smooth.quantile(rnbinom(100000, 2 * treat.SCL.linear.50.incorrect2incorrect + 1, expit(eval(parse(text = base.test)) +  treat.SCL.linear.50.incorrect2incorrect * (1 + x.test[, 2] + x.test[, 3]^2))), 0.50)
  v.SCL.gaussian.50.incorrect2incorrect <- smooth.quantile(rnbinom(100000, 2 * treat.SCL.gaussian.50.incorrect2incorrect + 1, expit(eval(parse(text = base.test)) +  treat.SCL.gaussian.50.incorrect2incorrect * (1 + x.test[, 2] + x.test[, 3]^2))), 0.50)

  return(list(
    SCL.linear.correctcorrect = c(v.SCL.linear.25.correctcorrect, mr.SCL.linear.25.correctcorrect, v.SCL.linear.50.correctcorrect, mr.SCL.linear.50.correctcorrect),
    SCL.gaussian.correctcorrect = c(v.SCL.gaussian.25.correctcorrect, mr.SCL.gaussian.25.correctcorrect, v.SCL.gaussian.50.correctcorrect, mr.SCL.gaussian.50.correctcorrect),
    SCL.linear.incorrect1correct = c(v.SCL.linear.25.incorrect1correct, mr.SCL.linear.25.incorrect1correct, v.SCL.linear.50.incorrect1correct, mr.SCL.linear.50.incorrect1correct),
    SCL.gaussian.incorrect1correct = c(v.SCL.gaussian.25.incorrect1correct, mr.SCL.gaussian.25.incorrect1correct, v.SCL.gaussian.50.incorrect1correct, mr.SCL.gaussian.50.incorrect1correct),
    SCL.linear.incorrect2correct = c(v.SCL.linear.25.incorrect2correct, mr.SCL.linear.25.incorrect2correct, v.SCL.linear.50.incorrect2correct, mr.SCL.linear.50.incorrect2correct),
    SCL.gaussian.incorrect2correct = c(v.SCL.gaussian.25.incorrect2correct, mr.SCL.gaussian.25.incorrect2correct, v.SCL.gaussian.50.incorrect2correct, mr.SCL.gaussian.50.incorrect2correct),
    SCL.linear.correctincorrect = c(v.SCL.linear.25.correctincorrect, mr.SCL.linear.25.correctincorrect, v.SCL.linear.50.correctincorrect, mr.SCL.linear.50.correctincorrect),
    SCL.gaussian.correctincorrect = c(v.SCL.gaussian.25.correctincorrect, mr.SCL.gaussian.25.correctincorrect, v.SCL.gaussian.50.correctincorrect, mr.SCL.gaussian.50.correctincorrect),
    SCL.linear.incorrect1incorrect = c(v.SCL.linear.25.incorrect1incorrect, mr.SCL.linear.25.incorrect1incorrect, v.SCL.linear.50.incorrect1incorrect, mr.SCL.linear.50.incorrect1incorrect),
    SCL.gaussian.incorrect1incorrect = c(v.SCL.gaussian.25.incorrect1incorrect, mr.SCL.gaussian.25.incorrect1incorrect, v.SCL.gaussian.50.incorrect1incorrect, mr.SCL.gaussian.50.incorrect1incorrect),
    SCL.linear.incorrect2incorrect = c(v.SCL.linear.25.incorrect2incorrect, mr.SCL.linear.25.incorrect2incorrect, v.SCL.linear.50.incorrect2incorrect, mr.SCL.linear.50.incorrect2incorrect),
    SCL.gaussian.incorrect2incorrect = c(v.SCL.gaussian.25.incorrect2incorrect, mr.SCL.gaussian.25.incorrect2incorrect, v.SCL.gaussian.50.incorrect2incorrect, mr.SCL.gaussian.50.incorrect2incorrect)))
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

# calculate the true optimal quantile 
case3.q.25 <- uniroot(opt.q.case3, c(-100, 100), tau = 0.25)$root
case3.q.50 <- uniroot(opt.q.case3, c(-100, 100), tau = 0.50)$root

# calculate the true optimal treatment regime for testing set
case3.opt.q.25 <- opt.reg.case3(case3.q.25)
case3.opt.q.50 <- opt.reg.case3(case3.q.50)

# parallel setup
cl <- makeCluster(getOption("cl.cores", clnum))
clusterEvalQ(cl, {
  library(MASS)
  library(glmnet)
  library(WeightSVM)
  library(quantoptr)
  library(quantreg)
})
clusterExport(cl, varlist = list("SCL.count", "smooth.quantile", "expit"), envir = .GlobalEnv)

# parallel computing
system.time({
  res <- parLapply(cl, 1:times, simulation, n = n, d=d, base = case3.base, base.test = case3.base.test, x.test = x, opt.q.25 = case3.opt.q.25, opt.q.50 = case3.opt.q.50)
})
stopCluster(cl)

# extract results
# value function and misclassification rates
case3.SCL.linear.correctcorrect <- array(0, dim = c(times, 4))
case3.SCL.gaussian.correctcorrect <- array(0, dim = c(times, 4))


case3.SCL.linear.incorrect1correct <- array(0, dim = c(times, 4))
case3.SCL.gaussian.incorrect1correct <- array(0, dim = c(times, 4))
case3.SCL.linear.incorrect2correct <- array(0, dim = c(times, 4))
case3.SCL.gaussian.incorrect2correct <- array(0, dim = c(times, 4))

case3.SCL.linear.correctincorrect <- array(0, dim = c(times, 4))
case3.SCL.gaussian.correctincorrect <- array(0, dim = c(times, 4))

case3.SCL.linear.incorrect1incorrect <- array(0, dim = c(times, 4))
case3.SCL.gaussian.incorrect1incorrect <- array(0, dim = c(times, 4))
case3.SCL.linear.incorrect2incorrect <- array(0, dim = c(times, 4))
case3.SCL.gaussian.incorrect2incorrect <- array(0, dim = c(times, 4))


for (j in 1:times) {
  case3.SCL.linear.correctcorrect[j, ] <- res[[j]]$SCL.linear.correctcorrect
  case3.SCL.gaussian.correctcorrect[j, ] <- res[[j]]$SCL.gaussian.correctcorrect

  case3.SCL.linear.incorrect1correct[j, ] <- res[[j]]$SCL.linear.incorrect1correct
  case3.SCL.gaussian.incorrect1correct[j, ] <- res[[j]]$SCL.gaussian.incorrect1correct
  case3.SCL.linear.incorrect2correct[j, ] <- res[[j]]$SCL.linear.incorrect2correct
  case3.SCL.gaussian.incorrect2correct[j, ] <- res[[j]]$SCL.gaussian.incorrect2correct

  case3.SCL.linear.correctincorrect[j, ] <- res[[j]]$SCL.linear.correctincorrect
  case3.SCL.gaussian.correctincorrect[j, ] <- res[[j]]$SCL.gaussian.correctincorrect

  case3.SCL.linear.incorrect1incorrect[j, ] <- res[[j]]$SCL.linear.incorrect1incorrect
  case3.SCL.gaussian.incorrect1incorrect[j, ] <- res[[j]]$SCL.gaussian.incorrect1incorrect
  case3.SCL.linear.incorrect2incorrect[j, ] <- res[[j]]$SCL.linear.incorrect2incorrect
  case3.SCL.gaussian.incorrect2incorrect[j, ] <- res[[j]]$SCL.gaussian.incorrect2incorrect
}
