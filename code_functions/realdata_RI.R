library(speff2trial)
library(MASS)
library(glmnet)
library(WeightSVM)
library(quantoptr)
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

fold=5
ri.results <- foreach(i = 1:times, .combine = "rbind", .packages = c("glmnet","MASS","WeightSVM", "quantoptr"), .export = c("predict.QIQlearningcontinuous", "predict.glmnetformula")) %dopar% { 
  #construct training and testing data
  set.seed(i)
  index <- sample(1:nrow(data_real), nrow(data_real), replace = FALSE)
  index1 <- index[1:(2*round(nrow(data_real)/fold))]
  index2 <- index[(2*round(nrow(data_real)/fold)+1):(4*round(nrow(data_real)/fold))]
  index3 <- index[(4*round(nrow(data_real)/fold)+1):nrow(data_real)]

  x_train1 <- x[index1,]
  y_train1 <- y[index1]
  a_train1 <- a[index1]
  data1 <- list(y=y_train1, x=x_train1, a=a_train1)

  x_train2 <- x[index2,]
  y_train2 <- y[index2]
  a_train2 <- a[index2]
  data2 <- list(y=y_train2, x=x_train2, a=a_train2)

  x.test <- x[index3,]
  y.test <- y[index3]
  a.test <- a[index3]

  ### estimate the quantile OTRs in the first training set for different methods
  ## estimation of nuisance functions
  # estimation of the survival model
  d = ncol(x_train1)-1
  dataframe <- data.frame(y = y_train1, x = unname(x_train1[,-1]), a = a_train1)
  formula <-"y ~ "
  for(i in 1:d){
    formula = paste0(formula, paste("x.", i, sep = ""))
    formula = paste0(formula, " + ")
    formula = paste0(formula, paste("I(a*x.", i, ")", sep = ""))
    formula = paste0(formula, " + ")
  }
  formula <- as.formula(paste0(formula, "a"))
  xa_train1 <- model.matrix(formula, dataframe)[, -1]

  fit.f.y.ax1 <- glmnet(xa_train1, y_train1, lambda=0)
  f.y.ax1 <- list(fit=fit.f.y.ax1, formula=formula)
  class(f.y.ax1) <- "glmnetformula"

  r <- y_train1 - predict(fit.f.y.ax1, xa_train1)
  fit.f.r.ax1 <- glmnet(xa_train1, r^2, family = Gamma(link = "log"), lambda=0)
  f.r.ax1 <- list(fit=fit.f.r.ax1, formula=formula)
  class(f.r.ax1) <- "glmnetformula"

  # the propensity score is fixed as 0.48 according to Wang 2018
  mean.a <- 0.48

  ### Successive classification learning
  data1 <- list(y=y_train1, x=unname(x_train1), a=a_train1)
  n <- length(y_train1)
  kappa <- sd(y_train1) / sqrt(n) / 6
  epsilon <- 1 / sqrt(n) / 2

  res.SCL.linear.25.1 <- SCL.continuous(data1, 0.25, kappa, epsilon, mean.a, f.y.ax1, f.r.ax1,kernel="linear")
  res.SCL.gaussian.25.1 <- SCL.continuous(data1,0.25, kappa, epsilon, mean.a, f.y.ax1, f.r.ax1, ,kernel="gaussian")

  res.SCL.linear.50.1 <- SCL.continuous(data1,0.50, kappa, epsilon, mean.a, f.y.ax1, f.r.ax1,kernel="linear")
  res.SCL.gaussian.50.1 <- SCL.continuous(data1,0.50, kappa, epsilon, mean.a, f.y.ax1, f.r.ax1,kernel="gaussian")

  res.SCL.linear.75.1 <- SCL.continuous(data1,0.75, kappa, epsilon, mean.a, f.y.ax1, f.r.ax1,kernel="linear")
  res.SCL.gaussian.75.1 <- SCL.continuous(data1,0.75, kappa, epsilon, mean.a, f.y.ax1, f.r.ax1,kernel="gaussian")

  # QI-learning method
  res.QIQ.25.1 <- QIQ.learning.continuous(data1, 0.25, f.y.ax1, f.r.ax1)
  res.QIQ.50.1 <- QIQ.learning.continuous(data1, 0.50, f.y.ax1, f.r.ax1)
  res.QIQ.75.1 <- QIQ.learning.continuous(data1, 0.75, f.y.ax1, f.r.ax1)

  # wang's method
  data.wang <- data.frame(y=y_train1, x=x_train1[,-1], a=a_train1)
  regimeClass <- as.formula(a ~ .-y)
  moCondQuant_0 <- as.formula(y ~ .-a)
  moCondQuant_1 <- as.formula(y ~ .-a)
  moPropen <- as.formula(a ~ 1)

  res.wang.50.1 <- DR_Qopt(data=data.wang, regimeClass = regimeClass, tau = 0.5,
                 moPropen = moPropen ,
                 moCondQuant_0 = moCondQuant_0,
                 moCondQuant_1 = moCondQuant_1, it.num = 2, pop.size = 500
                 )
  res.wang.25.1 <- DR_Qopt(data=data.wang, regimeClass = regimeClass, tau = 0.25,
                  moPropen = moPropen ,
                  moCondQuant_0 = moCondQuant_0,
                  moCondQuant_1 = moCondQuant_1, it.num = 2, pop.size = 500
                  )
  res.wang.75.1 <- DR_Qopt(data=data.wang, regimeClass = regimeClass, tau = 0.75,
                  moPropen = moPropen ,
                  moCondQuant_0 = moCondQuant_0,
                  moCondQuant_1 = moCondQuant_1, it.num = 2, pop.size = 500
                  )

  ### estimate the quantile OTRs in the second training set for different methods
  # estimation of the survival model
  dataframe <- data.frame(y = y_train2, x = unname(x_train2[,-1]), a = a_train2)
  formula <-"y ~ "
  for(i in 1:d){
    formula = paste0(formula, paste("x.", i, sep = ""))
    formula = paste0(formula, " + ")
    formula = paste0(formula, paste("I(a*x.", i, ")", sep = ""))
    formula = paste0(formula, " + ")
  }
  formula <- as.formula(paste0(formula, "a"))
  xa_train2 <- model.matrix(formula, dataframe)[, -1]

  fit.f.y.ax2 <- glmnet(xa_train2, y_train2, lambda=0)
  f.y.ax2 <- list(fit=fit.f.y.ax2, formula=formula)
  class(f.y.ax2) <- "glmnetformula"

  r <- y_train2 - predict(fit.f.y.ax2, xa_train2)
  fit.f.r.ax2 <- glmnet(xa_train2, r^2, family = Gamma(link = "log"), lambda=0)
  f.r.ax2 <- list(fit=fit.f.r.ax2, formula=formula)
  class(f.r.ax2) <- "glmnetformula"

  #successive classification learning
  data2 <- list(y=y_train2, x=unname(x_train2), a=a_train2)
  n2 <- length(y_train2)
  kappa <- sd(y_train2) / sqrt(n2) / 6
  epsilon <- 1 / sqrt(n2) / 2

  res.SCL.linear.25.2 <- SCL.continuous(data2,0.25, kappa, epsilon, mean.a, f.y.ax2, f.r.ax2, ,kernel="linear")
  res.SCL.gaussian.25.2 <- SCL.continuous(data2,0.25, kappa, epsilon, mean.a, f.y.ax2, f.r.ax2, ,kernel="gaussian")
  res.SCL.linear.50.2 <- SCL.continuous(data2,0.50, kappa, epsilon, mean.a, f.y.ax2, f.r.ax2, kernel="linear")
  res.SCL.gaussian.50.2 <- SCL.continuous(data2,0.50, kappa, epsilon, mean.a, f.y.ax2, f.r.ax2, kernel="gaussian")
  res.SCL.linear.75.2 <- SCL.continuous(data2,0.75, kappa, epsilon, mean.a, f.y.ax2, f.r.ax2, kernel="linear")
  res.SCL.gaussian.75.2 <- SCL.continuous(data2,0.75, kappa, epsilon, mean.a, f.y.ax2, f.r.ax2, kernel="gaussian")

  # QIQ learning
  res.QIQ.25.2 <- QIQ.learning.continuous(data2, 0.25, f.y.ax2, f.r.ax2)
  res.QIQ.50.2 <- QIQ.learning.continuous(data2, 0.50, f.y.ax2, f.r.ax2)
  res.QIQ.75.2 <- QIQ.learning.continuous(data2, 0.75, f.y.ax2, f.r.ax2)

  # wang's method
  data.wang <- data.frame(y=y_train2, x=x_train2[,-1], a=a_train2)
  regimeClass <- as.formula(a ~ .-y)
  moCondQuant_0 <- as.formula(y ~ .-a)
  moCondQuant_1 <- as.formula(y ~ .-a)
  moPropen <- as.formula(a ~ 1)

  res.wang.50.2 <- DR_Qopt(data=data.wang, regimeClass = regimeClass, tau = 0.5,
                  moPropen = moPropen ,
                 moCondQuant_0 = moCondQuant_0,
                 moCondQuant_1 = moCondQuant_1, it.num = 2, pop.size = 500
                 )
  res.wang.25.2 <- DR_Qopt(data=data.wang, regimeClass = regimeClass, tau = 0.25,
                  moPropen = moPropen ,
                  moCondQuant_0 = moCondQuant_0,
                  moCondQuant_1 = moCondQuant_1, it.num = 2, pop.size = 500
                  )
  res.wang.75.2 <- DR_Qopt(data=data.wang, regimeClass = regimeClass, tau = 0.75,
                  moPropen = moPropen ,
                  moCondQuant_0 = moCondQuant_0,
                  moCondQuant_1 = moCondQuant_1, it.num = 2, pop.size = 500
                  )
 
  # estimated quantile optimal treatment in test data
  x.test <- unname(x.test)
  treat.SCL.linear.25.1 <- as.numeric(predict(res.SCL.linear.25.1$opt.reg, data.frame(x=x.test[,-1])))-1
  treat.SCL.gaussian.25.1 <- as.numeric(predict(res.SCL.gaussian.25.1$opt.reg, data.frame(x=x.test[,-1])))-1
  treat.SCL.linear.50.1 <- as.numeric(predict(res.SCL.linear.50.1$opt.reg, data.frame(x=x.test[,-1])))-1
  treat.SCL.gaussian.50.1 <- as.numeric(predict(res.SCL.gaussian.50.1$opt.reg, data.frame(x=x.test[,-1])))-1
  treat.SCL.linear.75.1 <- as.numeric(predict(res.SCL.linear.75.1$opt.reg, data.frame(x=x.test[,-1])))-1
  treat.SCL.gaussian.75.1 <- as.numeric(predict(res.SCL.gaussian.75.1$opt.reg, data.frame(x=x.test[,-1])))-1

  treat.QIQ.25.1 <- predict(res.QIQ.25.1$opt.reg, x.test[,-1])
  treat.QIQ.50.1 <- predict(res.QIQ.50.1$opt.reg, x.test[,-1])
  treat.QIQ.75.1 <- predict(res.QIQ.75.1$opt.reg, x.test[,-1])

  if(sum(res.wang.25.1$coefficients)=='NaN'){
    treat.wang.25.1 <- x.test %*% rep(0,dim(x.test)[2])
  }else{
    treat.wang.25.1<- x.test %*% res.wang.25.1$coefficients>0
  }
  
  if(sum(res.wang.50.1$coefficients)=='NaN'){
    treat.wang.50.1 <- x.test %*% rep(0,dim(x.test)[2])
  }else{
    treat.wang.50.1<- x.test %*% res.wang.50.1$coefficients>0
  }

  if(sum(res.wang.75.1$coefficients)=='NaN'){
    treat.wang.75.1 <- x.test %*% rep(0,dim(x.test)[2])
  }else{
    treat.wang.75.1<- x.test %*% res.wang.75.1$coefficients>0
  }

  treat.SCL.linear.25.2 <- as.numeric(predict(res.SCL.linear.25.2$opt.reg, data.frame(x=x.test[,-1])))-1
  treat.SCL.gaussian.25.2 <- as.numeric(predict(res.SCL.gaussian.25.2$opt.reg, data.frame(x=x.test[,-1])))-1
  treat.SCL.linear.50.2 <- as.numeric(predict(res.SCL.linear.50.2$opt.reg, data.frame(x=x.test[,-1])))-1
  treat.SCL.gaussian.50.2 <- as.numeric(predict(res.SCL.gaussian.50.2$opt.reg, data.frame(x=x.test[,-1])))-1
  treat.SCL.linear.75.2 <- as.numeric(predict(res.SCL.linear.75.2$opt.reg, data.frame(x=x.test[,-1])))-1
  treat.SCL.gaussian.75.2 <- as.numeric(predict(res.SCL.gaussian.75.2$opt.reg, data.frame(x=x.test[,-1])))-1

  treat.QIQ.25.2 <- predict(res.QIQ.25.2$opt.reg, x.test[,-1])
  treat.QIQ.50.2 <- predict(res.QIQ.50.2$opt.reg, x.test[,-1])
  treat.QIQ.75.2 <- predict(res.QIQ.75.2$opt.reg, x.test[,-1])

  if(sum(res.wang.25.2$coefficients)=='NaN'){
    treat.wang.25.2 <- x.test %*% rep(0,dim(x.test)[2])
  }else{
    treat.wang.25.2<- x.test %*% res.wang.25.2$coefficients>0
  }
  
  if(sum(res.wang.50.2$coefficients)=='NaN'){
    treat.wang.50.2 <- x.test %*% rep(0,dim(x.test)[2])
  }else{
    treat.wang.50.2<- x.test %*% res.wang.50.2$coefficients>0
  }

  if(sum(res.wang.75.2$coefficients)=='NaN'){
    treat.wang.75.2 <- x.test %*% rep(0,dim(x.test)[2])
  }else{
    treat.wang.75.2<- x.test %*% res.wang.75.2$coefficients>0
  }

  # calculate the rand index for competing methods
  ri.SCL.linear.25 <- sum(treat.SCL.linear.25.1 == treat.SCL.linear.25.2)/length(y.test)
  ri.SCL.linear.50 <- sum(treat.SCL.linear.50.1 == treat.SCL.linear.50.2)/length(y.test)
  ri.SCL.linear.75 <- sum(treat.SCL.linear.75.1 == treat.SCL.linear.75.2)/length(y.test)

  ri.SCL.gaussian.25 <- sum(treat.SCL.gaussian.25.1 == treat.SCL.gaussian.25.2)/length(y.test)
  ri.SCL.gaussian.50 <- sum(treat.SCL.gaussian.50.1 == treat.SCL.gaussian.50.2)/length(y.test)
  ri.SCL.gaussian.75 <- sum(treat.SCL.gaussian.75.1 == treat.SCL.gaussian.75.2)/length(y.test)

  ri.QIQ.25 <- sum(treat.QIQ.25.1 == treat.QIQ.25.2)/length(y.test)
  ri.QIQ.50 <- sum(treat.QIQ.50.1 == treat.QIQ.50.2)/length(y.test)
  ri.QIQ.75 <- sum(treat.QIQ.75.1 == treat.QIQ.75.2)/length(y.test)

  ri.wang.25 <- sum(treat.wang.25.1 == treat.wang.25.2)/length(y.test)
  ri.wang.50 <- sum(treat.wang.50.1 == treat.wang.50.2)/length(y.test)
  ri.wang.75 <- sum(treat.wang.75.1 == treat.wang.75.2)/length(y.test)

  return(c("SCL.linear.25"=unname(ri.SCL.linear.25),
           "SCL.linear.50"=unname(ri.SCL.linear.50),
           "SCL.linear.75"=unname(ri.SCL.linear.75),
           "SCL.gaussian.25"=unname(ri.SCL.gaussian.25),
           "SCL.gaussian.50"=unname(ri.SCL.gaussian.50),
           "SCL.gaussian.75"=unname(ri.SCL.gaussian.75),
           "QIQ.25"=unname(ri.QIQ.25),
           "QIQ.50"=unname(ri.QIQ.50),
           "QIQ.75"=unname(ri.QIQ.75),
           "wang.25"=unname(ri.wang.25),
           "wang.50"=unname(ri.wang.50),
           "wang.75"=unname(ri.wang.75)
          )
         )
}

stopCluster(cl)

# extract the rand index 
realdata.ri.SCL.linear <- array(0, dim = c(times, 3))
realdata.ri.SCL.gaussian <- array(0, dim = c(times, 3))
realdata.ri.QIQ <- array(0, dim = c(times, 3))
realdata.ri.wang <- array(0, dim = c(times, 3))

for (j in 1:times) {
  realdata.ri.SCL.linear[j, ] <- ri.results[j, c("SCL.linear.25","SCL.linear.50","SCL.linear.75")]
  realdata.ri.SCL.gaussian[j, ] <- ri.results[j, c("SCL.gaussian.25","SCL.gaussian.50","SCL.gaussian.75")]
  realdata.ri.QIQ[j, ] <- ri.results[j, c("QIQ.25","QIQ.50","QIQ.75")]
  realdata.ri.wang[j, ] <- ri.results[j, c("wang.25","wang.50","wang.75")]
}