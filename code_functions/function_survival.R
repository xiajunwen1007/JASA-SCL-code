### Successive classification learning method for estimating quantile optimal treatment regimes for survival data.
#' @param data A list containing the observed outcome `tildey`, censoring indicator `delta`, treatment indicator `a`, and covariates `x`.
#' @param tau The quantile level for which the quantile optimal treatment regime is to be estimated.
#' @param kappa A small positive number to control the precision of the binary search procedure.
#' @param epsilon A small positive number to control the convergence criterion of the binary search procedure.
#' @param mean.a The propensity score pi(A=1|X).
#' @param f.y.x1 The conditional survival model of the survival time f.y.x1 for A=1.
#' @param f.y.x0 The conditional survival model of the survival time f.y.x0 for A=0.
#' @param f.c.ax The conditional survival model of the censoring time f.c.ax.
#' @param kernel The kernel type to construct quantile optimal treatment regimes, which can be "gaussian" and "linear". Default is "linear".
#' @param min The minimum value for the binary search procedure.
#' @param max The maximum value for the binary search procedure.
#' @param h The bandwidth for the soft regime Phi
#' @return A list containing opt.reg and q, representing the estimated quantile OTR and estimated quantile, respectively.
SCL.survival <- function(data, tau, kappa=0, epsilon=0, mean.a, f.y.x1, f.y.x0, f.c.ax,  kernel="linear", min=0, max=0, h=0) {
  # extract data components
  order <- order(data$tildey)
  tildey <- data$tildey[order]
  delta <- data$delta[order]
  a <- data$a[order]
  x <- data$x[order, ]
  n <- nrow(x)

  # initialize the binary search procedure
  if (min == 0 && max == 0) {
    min <- quantile(tildey, probs = 0.05)
    max <- quantile(tildey, probs = 0.95)
  }
  mid <- (max + min) / 2

  # parameters will be used in the binary search procedure
  # the smoothed regime function
  Phi <- function(x,h) pnorm(x/h)

  # the bandwidth for the smoothed regime function
  if (h==0){
     h <- 0.2/log(n)
  }
 
  # kappa
  if (kappa == 0) {
    kappa <- sd(tildey) / sqrt(n) / 6
  }
  # epsilon
  if (epsilon == 0) {
    epsilon <- 1 / sqrt(n) / 2
  }

  # the binary search procedure
  for (iter in 1:100) {
    newx <- data.frame(tildey=tildey, delta=delta,x=x[,-1])
    newx1 <- data.frame(tildey=tildey, delta=delta,a=1,x=x[,-1])
    newx0 <- data.frame(tildey=tildey, delta=delta,a=0,x=x[,-1])

    times <- tildey[tildey<=mid]

    # conditional survival function at mid
    S1 <- pmax(t(predict(f.y.x1, newx, times)),0.05)
    S0 <- pmax(t(predict(f.y.x0, newx, times)),0.05)
    C1 <- pmax(t(predict(f.c.ax, newx1, times)),0.05)
    C0 <- pmax(t(predict(f.c.ax, newx0, times)),0.05)

    dC1 <-cbind(1-C1[, 1], C1[, -ncol(C1)]- C1[, -1])
    dC0 <- cbind(1-C0[, 1],C0[, -ncol(C0)]- C0[, -1])
    dNC <- sapply(1:length(times), function(i) as.numeric(tildey == times[i] & delta == 0))
    IYt <- sapply(1:length(times), function(i) as.numeric(tildey >= times[i]))

    SC1 <- pmax(sapply(1:n, function(i) predict(f.c.ax, newx1[i,], tildey[i])), 0.05)
    SC0 <- pmax(sapply(1:n, function(i) predict(f.c.ax, newx0[i,], tildey[i])), 0.05)

    # efficient influence function to achieve the doubly robust property
    psi1 <- as.numeric(tildey > mid) * a *delta / (mean.a*SC1) + (1 - a / mean.a) * c(predict(f.y.x1, newx, mid))  + a* c(predict(f.y.x1, newx, mid))/mean.a*rowSums(dNC/C1/S1-IYt*dC1/(C1)^2/S1)  
    psi0 <-   as.numeric(tildey > mid) * (1 - a) * delta / ((1 - mean.a)*SC0) +(1-(1-a)/(1-mean.a)) * c(predict(f.y.x0, newx, mid)) + (1-a)* c(predict(f.y.x0, newx, mid))/(1-mean.a)*rowSums(dNC/C0/S0-IYt*dC0/(C0)^2/S0) 

    # calculate the new treatment indicator and the weights
    newa <- sign(psi1 - psi0)
    w <- abs(psi1 - psi0)
    data.regime <- data.frame(y = factor(newa), x = x[,-1])

    # if the kernel is "gaussian", use the gaussian kernel, otherwise use the linear kernel
    if (kernel == "gaussian") {
      # tune the cost parameter for the weighted SVM
       tn <- tune_wsvm(y ~ ., data = data.regime, kernel = "radial", ranges = list(cost = c(0.025 / 2, 0.05 / 2, 0.1 / 2, 1 / 2), gamma = c(0.025 / 2, 0.05 / 2, 0.1 / 2, 1 / 2)), weight = w)
    
      # fit the weighted SVM model
      opt.reg <- wsvm(y ~ ., data = data.regime, kernel = "radial", cost = tn$best.parameters$cost, gamma = tn$best.parameters$gamma, weight = w,probability=TRUE)

      # calculate the optimal treatment regime
      regime = Phi(newa[1]*attr(predict(opt.reg, data.regime[,-1], decision.values=TRUE),"decision.values"),h)
    } else {
      # use the linear kernel
      tn <- tune_wsvm(y ~ ., data = data.regime, kernel = "linear", ranges = list(cost = c( 0.05 / 2, 0.1 / 2, 1 / 2)), weight = w)
      opt.reg <- wsvm(y ~ ., data = data.regime, kernel = "linear", cost = tn$best.parameters$cost, weight = w)
      regime = Phi(newa[1]*attr(predict(opt.reg, data.regime[,-1], decision.values=TRUE),"decision.values"),h)
    }

    # stop the binary search procedure or update the search range
    if (max - min <= kappa || abs(mean(regime* (psi1-psi0) + psi0) - 1 + tau) <= epsilon) {
      break
    }
    if( mean( regime* psi1 + (1-regime)*psi0)  >= 1 - tau) {
      min <- mid 
    } else {
      max <- mid 
    }

    mid <- (max + min) / 2 
  }

  return(list(opt.reg=opt.reg, q=mid))
}

### functions to calculate the optimal treatment that maximizes the survival function at point q.
opt.reg.case1 <- function(q) {
  logq <- 4*log(q)
  return(as.numeric((logq - eval(parse(text = case1.base)) ) / 1 - ((logq - eval(parse(text = case1.base)) - x[, 2] -  x[, 3]) / 2) > 0))
}

### functions to calculate the difference between the optimal survival function at q and 1-tau. 
opt.q.case1 <- function(q, tau) {
  opt.a <- opt.reg.case1(q)
  logq <- 4*log(q)
  return(mean(pnorm(logq, mean = eval(parse(text = case1.base)) +  opt.a * (x[, 2] + x[, 3]), sd = opt.a+1, lower.tail = FALSE)) - 1 + tau)
}

### A prediction function for lognorm objects
### Input the lognorm object, new data, and time points to predict.
### Output the predicted values.
predict.lognorm <- function(fit, newdata, t) {
  X_new <- as.matrix(model.matrix(fit$terms, newdata))
  beta <- coef(fit)
  sigma <- fit$scale
  mu <- as.numeric(X_new %*% beta)
  return(1 - pnorm(outer(log(t), mu, "-") / sigma))
}
