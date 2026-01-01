# import necessary libraries
library(quantoptr)
library(doParallel)
library(foreach)

# Source the necessary functions
source("./code_functions/function_nonconvexity.R")

# Specify the data generation process. Note the first column of x is intercept.
base <- "-0.5-2*x[,2]+0.5*x[,3]+0.5*x[,4]+0.5*x[,5]+0.5*x[,6]"
base.test <- "-0.5-2*x.test[,2]+0.5*x.test[,3]+0.5*x.test[,4]+0.5*x.test[,5]+0.5*x.test[,6]"
err <- "rnorm(n)"
err.test <- "rnorm(m)"
n <- 500
d <- 5

# generate test data
set.seed(0)
m <- 100000 # sample size for test data
x.test <- cbind(1, array(rnorm(m * d), c(m, d)))
eps.test <- eval(parse(text = err.test))

# calculate the true optimal treatment regime at tau=0.5
q.50 <- uniroot(opt.q, c(-100, 100), tau = 0.50)$root
opt.q.50 <- opt.reg(q.50)

# set up for parallel computing
cl <- makeCluster(clnum)
registerDoParallel(cl)

res <- foreach(l = 1:times, .combine = "rbind", .packages = c("quantoptr")) %dopar%{
  # set seed
  set.seed(l)

  # generate data
  x <- cbind(1, array(rnorm(n * d), c(n,d)))
  ps <- 0.5
  a <- rbinom(n, 1, ps)
  eps <- eval(parse(text = err))
  y <- eval(parse(text = base)) + a * (x[, 2]+x[, 3]+x[,4]+x[,5]+x[,6]) + (2 * a + 1) * eps

  ### grid search method
  data <- list(x = x, a = a, y = y)
  res.gridsearch.50 <- grid_serach(data, 0.50)

  # wang's method
  data.wang <- data.frame(y = y, x = x[, -1], a = a)
  regimeClass <- as.formula("a ~ x.1+x.2+x.3+x.4+x.5")
  moPropen <- as.formula(a~1)

  res.wang.50 <- IPWE_Qopt(
    data = data.wang, regimeClass = regimeClass, tau = 0.50,
    moPropen = moPropen,
    pop.size = 500, it.num = 2, p_level =0
  )
  
  # evaluate whether wangs method find the global optimal
  if((res.gridsearch.50$hatQ-res.wang.50$hatQ)<=0.00001){
    flag <- 1
  }else{
    flag <- 0
  }

  # estimated quantile optimal treatment for test data
  treat.wang.50 <- as.numeric(x.test %*% res.wang.50$coefficients > 0)
  treat.gridsearch.50 <- as.numeric(x.test %*% res.gridsearch.50$coef > 0)

  # calculate the estimated quantile for competing methods
  v.wang.50 <- quantile(eval(parse(text = base.test)) + treat.wang.50 * (x.test[, 2] +x.test[,3]+x.test[,4]+x.test[,5]+x.test[,6]) + (2 * treat.wang.50 + 1) * eps.test, 0.50)
  v.gridsearch.50 <- quantile(eval(parse(text = base.test)) + treat.gridsearch.50 * (x.test[, 2]+x.test[,3]+x.test[,4]+x.test[,5]+x.test[,6]) + (2 * treat.gridsearch.50 + 1) * eps.test, 0.50)

  MR.wang.50 <- 1 - mean(treat.wang.50==opt.q.50)
  MR.gridsearch.50 <- 1 - mean(treat.gridsearch.50==opt.q.50)

  return(c("flag"=unname(flag),
           "v.wang.50"=unname(v.wang.50),
            "v.gridsearch.50"=unname(v.gridsearch.50),
            "MR.wang.50"=unname(MR.wang.50),
            "MR.gridsearch.50"=unname(MR.gridsearch.50)))
}

flag <- rep(0,times)
wang.50 <- array(0, c(times, 2))
gridsearch.50 <- array(0, c(times, 2))

for(i in 1:times){
  flag[i] <- res[i,"flag"]
  wang.50[i,1] <- res[i,"v.wang.50"]
  wang.50[i,2] <- res[i,"MR.wang.50"]
  gridsearch.50[i,1] <- res[i,"v.gridsearch.50"]
  gridsearch.50[i,2] <- res[i,"MR.gridsearch.50"]
}