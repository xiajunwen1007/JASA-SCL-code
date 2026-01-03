library(speff2trial)
library(MASS)
library(glmnet)
library(WeightSVM)
library(quantoptr)
library(quantreg)
library(doParallel)
library(foreach)

### source necessary functions
source("./code_functions/function_main.R")

### data preparation
data_real=rbind(ACTG175[ACTG175$arms==1,], ACTG175[ACTG175$arms==3,])
a=data_real$arms
for(i in 1:length(a)){
  if(a[i]==3){
    a[i]=0
  }else{
    a[i]=1
  }
}
y=data_real$cd420 
x=cbind(1,scale(as.matrix(data_real[,c(-1,-8,-9,-10,-11,-15,-17,-18,-20,-21,-22,-24,-25,-26,-27)])))

### setup for parallel computing
cl <- makeCluster(clnum)
registerDoParallel(cl)

fold <- 5
value.results <- foreach(i = 1:times, .combine = "rbind", .packages = c("glmnet","MASS","WeightSVM", "quantreg", "quantoptr"), .export=c("predict.QIQlearningcontinuous", "predict.glmnetformula")) %dopar% { 
  set.seed(i)
  # construct training and testing data
  index <- sample(1:nrow(data_real), nrow(data_real), replace = FALSE)
  index1 <- index[1:(4*round(nrow(data_real)/fold))]
  index2 <- index[(4*round(nrow(data_real)/fold)+1):nrow(data_real)]

  x_train <- x[index2,]
  y_train <- y[index2]
  a_train <- a[index2]
  x.test <- x[index1,]
  y.test <- y[index1]
  a.test <- a[index1]

  ### estimation of nuisance functions
  # estimation of the survival model
  d = ncol(x_train)-1
  dataframe <- data.frame(y = y_train, x = unname(x_train[,-1]), a = a_train)
  formula <-"y ~ "
  for(i in 1:d){
    formula = paste0(formula, paste("x.", i, sep = ""))
    formula = paste0(formula, " + ")
    formula = paste0(formula, paste("I(a*x.", i, ")", sep = ""))
    formula = paste0(formula, " + ")
  }
  formula <- as.formula(paste0(formula, "a"))
  xa_train <- model.matrix(formula, dataframe)[, -1]

  fit.f.y.ax <- glmnet(xa_train, y_train, lambda=0)
  f.y.ax <- list(fit=fit.f.y.ax, formula=formula)
  class(f.y.ax) <- "glmnetformula"

  r <- y_train - predict(fit.f.y.ax, xa_train)
  fit.f.r.ax <- glmnet(xa_train, r^2, family = Gamma(link = "log"), lambda=0)
  f.r.ax <- list(fit=fit.f.r.ax, formula=formula)
  class(f.r.ax) <- "glmnetformula"

  # the propensity score is fixed as 0.48 according to Wang 2018
  mean.a <- 0.48

  ### Successive classification learning
  data <- list(y=y_train, x=unname(x_train), a=a_train)
  n <- length(y_train)
  kappa <- sd(y_train) / sqrt(n) / 6
  epsilon <- 1 / sqrt(n) / 2

  res.SCL.linear.25 <- SCL.continuous(data,0.25, kappa, epsilon, mean.a, f.y.ax, f.r.ax, kernel="linear")
  res.SCL.gaussian.25 <- SCL.continuous(data,0.25, kappa, epsilon, mean.a, f.y.ax, f.r.ax, kernel="gaussian")
  res.SCL.linear.50 <- SCL.continuous(data,0.50, kappa, epsilon, mean.a, f.y.ax, f.r.ax, kernel="linear")
  res.SCL.gaussian.50 <- SCL.continuous(data,0.50, kappa, epsilon, mean.a, f.y.ax, f.r.ax, kernel="gaussian")
  res.SCL.linear.75 <- SCL.continuous(data,0.75, kappa, epsilon, mean.a, f.y.ax, f.r.ax, kernel="linear")
  res.SCL.gaussian.75 <- SCL.continuous(data,0.75, kappa, epsilon, mean.a, f.y.ax, f.r.ax, kernel="gaussian")

  # QIQ learning model
  res.QIQ.25 <- QIQ.learning.continuous(data, 0.25, f.y.ax, f.r.ax)
  res.QIQ.50 <- QIQ.learning.continuous(data, 0.5, f.y.ax, f.r.ax)
  res.QIQ.75 <- QIQ.learning.continuous(data, 0.75, f.y.ax, f.r.ax)

  # Wang's method
  data.wang <- data.frame(y=y_train, x=x_train[,-1], a=a_train)
  regimeClass <- as.formula(a ~ .-y)
  moCondQuant_0 <- as.formula(y ~ .-a)
  moCondQuant_1 <- as.formula(y ~ .-a)
  moPropen <- as.formula(a ~ 1)

  res.wang.25 <- DR_Qopt(data=data.wang, regimeClass = regimeClass, tau = 0.25,
                  moPropen = moPropen ,
                  moCondQuant_0 = moCondQuant_0,
                  moCondQuant_1 = moCondQuant_1, it.num = 2, pop.size = 500
                  )
  res.wang.50 <- DR_Qopt(data=data.wang, regimeClass = regimeClass, tau = 0.5,
                  moPropen = moPropen ,
                  moCondQuant_0 = moCondQuant_0,
                  moCondQuant_1 = moCondQuant_1, it.num = 2, pop.size = 500
                  )
  res.wang.75 <- DR_Qopt(data=data.wang, regimeClass = regimeClass, tau = 0.75,
                  moPropen = moPropen ,
                  moCondQuant_0 = moCondQuant_0,
                  moCondQuant_1 = moCondQuant_1, it.num = 2, pop.size = 500
                  )

  # estimated quantile optimal treatment for test data
  x.test <- unname(x.test)
  treat.SCL.linear.25 <- as.numeric(predict(res.SCL.linear.25$opt.reg, data.frame(x=x.test[,-1])))-1
  treat.SCL.gaussian.25 <- as.numeric(predict(res.SCL.gaussian.25$opt.reg, data.frame(x=x.test[,-1])))-1
  treat.SCL.linear.50 <- as.numeric(predict(res.SCL.linear.50$opt.reg, data.frame(x=x.test[,-1])))-1
  treat.SCL.gaussian.50 <- as.numeric(predict(res.SCL.gaussian.50$opt.reg, data.frame(x=x.test[,-1])))-1
  treat.SCL.linear.75 <- as.numeric(predict(res.SCL.linear.75$opt.reg, data.frame(x=x.test[,-1])))-1
  treat.SCL.gaussian.75 <- as.numeric(predict(res.SCL.gaussian.75$opt.reg, data.frame(x=x.test[,-1])))-1

  treat.QIQ.25 <- predict(res.QIQ.25$opt.reg, x.test[,-1])
  treat.QIQ.50 <- predict(res.QIQ.50$opt.reg, x.test[,-1])
  treat.QIQ.75 <- predict(res.QIQ.75$opt.reg, x.test[,-1])

  if(sum(res.wang.25$coefficients)=='NaN'){
    treat.wang.25 <- x.test %*% rep(0,dim(x.test)[2])
  }else{
    treat.wang.25<- x.test %*% res.wang.25$coefficients>0
  }

  if(sum(res.wang.50$coefficients)=='NaN'){
    treat.wang.50 <- x.test %*% rep(0,dim(x.test)[2])
  }else{
    treat.wang.50<- x.test %*% res.wang.50$coefficients>0
  }

  if(sum(res.wang.75$coefficients)=='NaN'){
    treat.wang.75 <- x.test %*% rep(0,dim(x.test)[2])
  }else{
    treat.wang.75<- x.test %*% res.wang.75$coefficients>0
  }

  # calculate the estimated quantile for competing methods
  v.SCL.linear.25 <- rq(y.test ~ 1, tau = 0.25, weights=as.numeric(treat.SCL.linear.25==a.test)/(a.test*mean.a+(1-a.test)*(1-mean.a)))$coefficients
  v.SCL.gaussian.25 <- rq(y.test ~ 1, tau = 0.25, weights=as.numeric(treat.SCL.gaussian.25==a.test)/(a.test*mean.a+(1-a.test)*(1-mean.a)))$coefficients
  v.SCL.linear.50 <- rq(y.test ~ 1, tau = 0.50, weights=as.numeric(treat.SCL.linear.50==a.test)/(a.test*mean.a+(1-a.test)*(1-mean.a)))$coefficients
  v.SCL.gaussian.50 <- rq(y.test ~ 1, tau = 0.50, weights=as.numeric(treat.SCL.gaussian.50==a.test)/(a.test*mean.a+(1-a.test)*(1-mean.a)))$coefficients
  v.SCL.linear.75 <- rq(y.test ~ 1, tau = 0.75, weights=as.numeric(treat.SCL.linear.75==a.test)/(a.test*mean.a+(1-a.test)*(1-mean.a)))$coefficients
  v.SCL.gaussian.75 <- rq(y.test ~ 1, tau = 0.75, weights=as.numeric(treat.SCL.gaussian.75==a.test)/(a.test*mean.a+(1-a.test)*(1-mean.a)))$coefficients

  v.QIQ.25 <- rq(y.test ~ 1, tau = 0.25, weights=as.numeric(treat.QIQ.25==a.test)/(a.test*mean.a+(1-a.test)*(1-mean.a)))$coefficients
  v.QIQ.50 <- rq(y.test ~ 1, tau = 0.50, weights=as.numeric(treat.QIQ.50==a.test)/(a.test*mean.a+(1-a.test)*(1-mean.a)))$coefficients
  v.QIQ.75 <- rq(y.test ~ 1, tau = 0.75, weights=as.numeric(treat.QIQ.75==a.test)/(a.test*mean.a+(1-a.test)*(1-mean.a)))$coefficients

  v.wang.25 <- rq(y.test ~ 1, tau = 0.25, weights=as.numeric(treat.wang.25==a.test)/(a.test*mean.a+(1-a.test)*(1-mean.a)))$coefficients
  v.wang.50 <- rq(y.test ~ 1, tau = 0.50, weights=as.numeric(treat.wang.50==a.test)/(a.test*mean.a+(1-a.test)*(1-mean.a)))$coefficients
  v.wang.75 <- rq(y.test ~ 1, tau = 0.75, weights=as.numeric(treat.wang.75==a.test)/(a.test*mean.a+(1-a.test)*(1-mean.a)))$coefficients

  return(c("SCL.linear.25"=unname(v.SCL.linear.25),
           "SCL.gaussian.25"=unname(v.SCL.gaussian.25),
           "QIQ.25"=unname(v.QIQ.25),
           "wang.25"=unname(v.wang.25),
           "SCL.linear.50"=unname(v.SCL.linear.50),
           "SCL.gaussian.50"=unname(v.SCL.gaussian.50),
           "QIQ.50"=unname(v.QIQ.50),
           "wang.50"=unname(v.wang.50),
           "SCL.linear.75"=unname(v.SCL.linear.75),
           "SCL.gaussian.75"=unname(v.SCL.gaussian.75),
           "QIQ.75"=unname(v.QIQ.75),
           "wang.75"=unname(v.wang.75)
         ))
}

stopCluster(cl)

# extract the value function 
realdata.value.SCL.linear <- array(0, dim = c(times, 3))
realdata.value.SCL.gaussian <- array(0, dim = c(times, 3))
realdata.value.QIQ <- array(0, dim = c(times, 3))
realdata.value.wang <- array(0, dim = c(times, 3))

for (j in 1:times) {
  realdata.value.SCL.linear[j, ] <- value.results[j, c("SCL.linear.25","SCL.linear.50","SCL.linear.75")]
  realdata.value.SCL.gaussian[j, ] <- value.results[j, c("SCL.gaussian.25","SCL.gaussian.50","SCL.gaussian.75")]
  realdata.value.QIQ[j, ] <- value.results[j, c("QIQ.25","QIQ.50","QIQ.75")]
  realdata.value.wang[j, ] <- value.results[j, c("wang.25","wang.50","wang.75")]
}