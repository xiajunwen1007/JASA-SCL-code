### import necessary libraries
library(MASS)
library(parallel)
library(WeightSVM)
library(quantoptr)

### load necessary functions
source("./code_functions/function_inconsistency.R")

### The main simulation function that generates data, fits different methods, and evaluates their performance.
### The input l denotes the index of the simulation run, the input n is the training sample size, and the input x.test is the test samples used to evaluate the methods. 
### The outputs include the value functions for different methods.
simulation <- function(l, n, x.test) {
  # generate data
  n <- n
  d <- 2

  x <- cbind(1, array(rnorm(n * d), c(n, d)))
  ps <- 0.5
  a <- rbinom(n, 1, ps)
  mean.y <- expit(-0.5+0.5 * x[, 2] + 0.5 * x[, 3]+a*(1.5+1.5*x[,2]-1.5*x[,3]))
  y <- rbinom(n, size = 1, prob = mean.y)

  data <- list(y = y, x = x, a = a)

  ### estimation of nuisance functions
  # estimation of survival function
  dataframe <- data.frame(y = y, x = x[, -1], a = a)
  formula <- as.formula("y ~ x.1 + x.2 + a + I(a * x.1) + I(a * x.2)")
  f.y.ax <- glm(formula = formula, data = dataframe, family = binomial(link = "logit"))

  # the propensity score is known as 0.50
  mean.a <- rep(0.5,n)

  # Successive classification learning
  kappa <- sd(y) / sqrt(n) / 6
  epsilon <- 1 / sqrt(n) / 2

  res.SCL.linear.40 <- SCL.binary(data, 0.40, kappa, epsilon, mean.a, f.y.ax, kernel="linear")
  res.SCL.linear.50 <- SCL.binary(data, 0.50, kappa, epsilon, mean.a, f.y.ax, kernel="linear")
  res.SCL.linear.60 <- SCL.binary(data, 0.60, kappa, epsilon, mean.a, f.y.ax, kernel="linear")

  # wang's method
  xformula <- ""
  for(j in 1:d){
    xformula <- paste(xformula,paste("x.", j, sep = ""),seq="")
    if (j < d) {
      xformula <- paste(xformula,"+",sep="")
    }
  }
  data.wang <- data.frame(y = y, x = x[, -1], a = a)
  regimeClass <- as.formula(paste("a ~", xformula))
  moPropen <- as.formula(a~1)

  res.wang.40 <- IPWE_Qopt(
    data = data.wang, regimeClass = regimeClass, tau = 0.40,
    moPropen = moPropen,
    pop.size = 500, it.num = 2
  )
  res.wang.50 <- IPWE_Qopt(
    data = data.wang, regimeClass = regimeClass, tau = 0.50,
    moPropen = moPropen,
    pop.size = 500, it.num = 2
  )
  res.wang.60 <- IPWE_Qopt(
    data = data.wang, regimeClass = regimeClass, tau = 0.60,
    moPropen = moPropen,
    pop.size = 500, it.num = 2
  )
  
  # evaluate the estimated treatment regimes on the test data
  treat.SCL.linear.40 <- as.numeric(predict(res.SCL.linear.40$opt.reg, x.test[, -1])) - 1
  treat.SCL.linear.50 <- as.numeric(predict(res.SCL.linear.50$opt.reg, x.test[, -1])) - 1
  treat.SCL.linear.60 <- as.numeric(predict(res.SCL.linear.60$opt.reg, x.test[, -1])) - 1

  if (sum(res.wang.40$coefficients) == "NaN") {
    treat.wang.40 <- x.test %*% rep(0, p)
  } else {
    treat.wang.40 <- x.test %*% res.wang.40$coefficients > 0
  }

  if (sum(res.wang.50$coefficients) == "NaN") {
    treat.wang.50 <- x.test %*% rep(0, p)
  } else {
    treat.wang.50 <- x.test %*% res.wang.50$coefficients > 0
  }

  if (sum(res.wang.60$coefficients) == "NaN") {
    treat.wang.60 <- x.test %*% rep(0, p)
  } else {
    treat.wang.60 <- x.test %*% res.wang.60$coefficients > 0
  }
  
  # calculate the value function for competing methods
  v.SCL.linear.40 <- as.numeric(mean(expit(-0.5+0.5 * x.test[, 2] + 0.5 * x.test[, 3]+ treat.SCL.linear.40*(1.5+1.5*x.test[,2]-1.5*x.test[,3]))-1+ 0.40)>0)
  v.SCL.linear.50 <- as.numeric(mean(expit(-0.5+0.5 * x.test[, 2] + 0.5 * x.test[, 3]+ treat.SCL.linear.50*(1.5+1.5*x.test[,2]-1.5*x.test[,3]))-1+ 0.50)>0)
  v.SCL.linear.60 <- as.numeric(mean(expit(-0.5+0.5 * x.test[, 2] + 0.5 * x.test[, 3]+ treat.SCL.linear.60*(1.5+1.5*x.test[,2]-1.5*x.test[,3]))-1+ 0.60)>0)

  v.wang.40 <- as.numeric(mean(expit(-0.5+0.5 * x.test[, 2] + 0.5 * x.test[, 3]+ treat.wang.40*(1.5+1.5*x.test[,2]-1.5*x.test[,3]))-1+ 0.40)>0)
  v.wang.50 <- as.numeric(mean(expit(-0.5+0.5 * x.test[, 2] + 0.5 * x.test[, 3]+ treat.wang.50*(1.5+1.5*x.test[,2]-1.5*x.test[,3]))-1+ 0.50)>0)
  v.wang.60 <- as.numeric(mean(expit(-0.5+0.5 * x.test[, 2] + 0.5 * x.test[, 3]+ treat.wang.60*(1.5+1.5*x.test[,2]-1.5*x.test[,3]))-1+ 0.60)>0)

  # calculate hatQ{Y(hat{d})}
  hatQ.wang.40 <- res.wang.40$hatQ
  hatQ.wang.50 <- res.wang.50$hatQ
  hatQ.wang.60 <- res.wang.60$hatQ

  return(list(
    v.40 = c(v.SCL.linear.40, v.wang.40),
    v.50 = c(v.SCL.linear.50, v.wang.50),
    v.60 = c(v.SCL.linear.60, v.wang.60),
    hatQ.wang.40 = hatQ.wang.40,
    hatQ.wang.50 = hatQ.wang.50,
    hatQ.wang.60 = hatQ.wang.60
  ))
}

### generate test data
d = 2
m = 100000 # sample size for test data
x <- cbind(1, array(rnorm(m * d), c(m,d))) # generate test data

### set up for parallel computing
cl <- makeCluster(getOption("cl.cores", clnum))
clusterEvalQ(cl, {
  library(MASS)
  library(glmnet)
  library(WeightSVM)
  library(quantoptr)
})
clusterExport(cl, varlist = list("expit", "SCL.binary"), envir = .GlobalEnv)

# parallel computing for sample size n=30,100,500,1000
clusterSetRNGStream(cl, 111)
system.time({
  res30 <- parLapply(cl, 1:times, simulation, n = 30, x.test = x)
})

clusterSetRNGStream(cl, 111)
system.time({
  res100 <- parLapply(cl, 1:times, simulation, n = 100, x.test = x)
})

clusterSetRNGStream(cl, 111)
system.time({
  res500 <- parLapply(cl, 1:times, simulation, n = 500, x.test = x)
})

clusterSetRNGStream(cl, 111)
system.time({
  res1000 <- parLapply(cl, 1:times, simulation, n = 1000, x.test = x)
})

stopCluster(cl)

# extract results
# extract value functions
v.SCL.linear.40 <- array(0, c(times, 4))
v.SCL.linear.50 <- array(0, c(times, 4))
v.SCL.linear.60 <- array(0, c(times, 4))

v.wang.40 <- array(0, c(times, 4))
v.wang.50 <- array(0, c(times, 4))
v.wang.60 <- array(0, c(times, 4))

# extract hatQ{Y(hat{d})}
hatQ.wang.40 <- array(0, c(times, 4))
hatQ.wang.50 <- array(0, c(times, 4))
hatQ.wang.60 <- array(0, c(times, 4))

for (i in 1:times) {
  v.SCL.linear.40[i, 1] <- res30[[i]]$v.40[1]
  v.SCL.linear.50[i, 1] <- res30[[i]]$v.50[1]
  v.SCL.linear.60[i, 1] <- res30[[i]]$v.60[1]

  v.SCL.linear.40[i, 2] <- res100[[i]]$v.40[1]
  v.SCL.linear.50[i, 2] <- res100[[i]]$v.50[1]
  v.SCL.linear.60[i, 2] <- res100[[i]]$v.60[1]

  v.SCL.linear.40[i, 3] <- res500[[i]]$v.40[1]
  v.SCL.linear.50[i, 3] <- res500[[i]]$v.50[1]
  v.SCL.linear.60[i, 3] <- res500[[i]]$v.60[1]

  v.SCL.linear.40[i, 4] <- res1000[[i]]$v.40[1]
  v.SCL.linear.50[i, 4] <- res1000[[i]]$v.50[1]
  v.SCL.linear.60[i, 4] <- res1000[[i]]$v.60[1]

  v.wang.40[i, 1] <- res30[[i]]$v.40[2]
  v.wang.50[i, 1] <- res30[[i]]$v.50[2]
  v.wang.60[i, 1] <- res30[[i]]$v.60[2]

  v.wang.40[i, 2] <- res100[[i]]$v.40[2]
  v.wang.50[i, 2] <- res100[[i]]$v.50[2]
  v.wang.60[i, 2] <- res100[[i]]$v.60[2]

  v.wang.40[i, 3] <- res500[[i]]$v.40[2]
  v.wang.50[i, 3] <- res500[[i]]$v.50[2]
  v.wang.60[i, 3] <- res500[[i]]$v.60[2]

  v.wang.40[i, 4] <- res1000[[i]]$v.40[2]
  v.wang.50[i, 4] <- res1000[[i]]$v.50[2]
  v.wang.60[i, 4] <- res1000[[i]]$v.60[2]

  hatQ.wang.40[i, 1] <- res30[[i]]$hatQ.wang.40
  hatQ.wang.50[i, 1] <- res30[[i]]$hatQ.wang.50
  hatQ.wang.60[i, 1] <- res30[[i]]$hatQ.wang.60

  hatQ.wang.40[i, 2] <- res100[[i]]$hatQ.wang.40
  hatQ.wang.50[i, 2] <- res100[[i]]$hatQ.wang.50
  hatQ.wang.60[i, 2] <- res100[[i]]$hatQ.wang.60

  hatQ.wang.40[i, 3] <- res500[[i]]$hatQ.wang.40
  hatQ.wang.50[i, 3] <- res500[[i]]$hatQ.wang.50
  hatQ.wang.60[i, 3] <- res500[[i]]$hatQ.wang.60

  hatQ.wang.40[i, 4] <- res1000[[i]]$hatQ.wang.40
  hatQ.wang.50[i, 4] <- res1000[[i]]$hatQ.wang.50
  hatQ.wang.60[i, 4] <- res1000[[i]]$hatQ.wang.60
}