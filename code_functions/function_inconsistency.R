### Successive classification learning method for estimating quantile optimal treatment regimes when outcomes are binary data.
#' @param data A list containing the outcome `y`, treatment indicator `a`, and covariates `x`.
#' @param tau The quantile level for which the quantile optimal treatment regime is to be estimated.
#' @param kappa A small positive number to control the precision of the binary search procedure.
#' @param epsilon A small positive number to control the convergence criterion of the binary search procedure.
#' @param mean.a The estimated propensity score pi(A=1|X).
#' @param f.y.ax The conditional survival model for a binary y.
#' @param kernel The kernel type to construct quantile optimal treatment regimes, which can be "gaussian" and "linear". Default is "linear".
#' @param min The minimum value for the binary search procedure.
#' @param max The maximum value for the binary search procedure.
#' @param h The bandwidth for the soft regime Phi
#' @return A list containing opt.reg and q, representing the estimated quantile OTR and estimated quantile, respectively.
SCL.binary <- function(data, tau, kappa=0, epsilon=0, mean.a, f.y.ax,  kernel="linear", min=0, max=0, h=0) {
  # extract data components
  y <- data$y
  a <- data$a
  x <- data$x
  n <- nrow(x)

  # the smoothed regime function
  Phi <-function(x,h) pnorm(x/h)
  # the bandwidth for the smoothed regime function
  h <- 0.2/log(n)
  # kappa
  if (kappa == 0) {
    kappa <- sd(y) / sqrt(n) / 6
  }

  # epsilon
  if (epsilon == 0) {
    epsilon <- 1 / sqrt(n) / 2
  }

  #initialize the binary search procedure
  if (min == 0 && max == 0) {
    min <- min(y)-1
    max <- max(y)
  }
  mid <- (max + min) / 2

  for (iter in 1:100) {
    newx1 <- data.frame(y=y,x=x[,-1], a=1)
    newx0 <- data.frame(y=y,x=x[,-1], a=0)

    # calculate the conditional smoothed survival function at mid
    qfloor <- floor(mid)
    qceiling <- ceiling(mid)
    if (qfloor == qceiling) {
      lambda <- 0
    } else {
      lambda <- (mid - qceiling) / (qfloor - qceiling)
    }

    # conditional smoothed survival function at mid
    S1 <- pbinom(qceiling, size=1, prob = predict(f.y.ax, newx1, type = "response"), lower.tail = FALSE) + lambda * dbinom(qceiling, size=1, prob = predict(f.y.ax, newx1, type = "response"))
    S0 <- pbinom(qceiling, size=1, prob = predict(f.y.ax, newx0, type = "response"), lower.tail = FALSE) + lambda * dbinom(qceiling, size=1, prob = predict(f.y.ax, newx0, type = "response"))
    
    # smoothed psi function values to achieve the doubly robust property
    psi1 <- S1 + (as.numeric(y > qceiling) + lambda * (as.numeric(y == qceiling)) - S1) / mean.a * a
    psi0 <- S0 + (as.numeric(y > qceiling) + lambda * (as.numeric(y == qceiling)) - S0) / (1 - mean.a) * (1 - a)

    # calculate the new treatment indicator and the weights
    newa <- sign(psi1 - psi0)
    w <- abs(psi1 - psi0)
    data.regime <- data.frame(y = factor(newa), x = x[,-1])
    
    # fit the weighted SVM model to obtain the estimated optimal treatment regime
    # When newa is all 1 and 0, create a model predict all x to newa[1]
    if(length(unique(newa))==1){
      regime <- newa[1]
      # Create a dummy data frame to train the model
      df <- data.frame(x = rnorm(10), y = rep(newa[1], 10))  # y is always newa[1]

      # Fit a regression model with no x-dependence
      opt.reg <- lm(y ~ 1, data = df) 
    }else{
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
        tn <- tune_wsvm(y ~ ., data = data.regime, kernel = "linear", ranges = list(cost = c(0.025 / 2, 0.05 / 2, 0.1 / 2, 1 / 2)), weight = w)
        opt.reg <- wsvm(y ~ ., data = data.regime, kernel = "linear", cost = tn$best.parameters$cost, weight = w)
        regime = Phi(newa[1]*attr(predict(opt.reg, data.regime[,-1], decision.values=TRUE),"decision.values"),h)
      }
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

expit <- function(x) {
  return(1 / (1 + exp(-x)))
}

### functions to calculate the difference between the optimal survival function at q and 1-tau.
opt.q <- function(q, tau) {
  opt.a <- opt.reg(q)

  qfloor <- floor(q)
  qceiling <- ceiling(q)

  if (qfloor == qceiling) {
    lambda <- 0
  } else {
    lambda <- (q - qceiling) / (qfloor - qceiling)
  }
  return(mean(pbinom(qceiling, size=1, prob = expit(-0.5+0.5 * x[, 2] + 0.5 * x[, 3]+opt.a*(1.5+1.5*x[,2]-1.5*x[,3])), lower.tail = FALSE) + lambda * dbinom(qceiling, size=1, prob =expit(-0.5+0.5 * x[, 2] + 0.5 * x[, 3]+opt.a*(1.5+1.5*x[,2]-1.5*x[,3])))) -1 + tau)
}

### functions to calculate the optimal treatment that maximizes the survival function at point q.
opt.reg <- function(q) {
  qfloor <- floor(q)
  qceiling <- ceiling(q)

  if (qfloor == qceiling) {
    lambda <- 0
  } else {
    lambda <- (q - qceiling) / (qfloor - qceiling)
  }

  S1 <- pbinom(qceiling, size=1, prob = expit(-0.5+0.5 * x[, 2] + 0.5 * x[, 3]+1*(1.5+1.5*x[,2]-1.5*x[,3])), lower.tail = FALSE) + lambda * dbinom(qceiling, size=1, prob =expit(-0.5+0.5 * x[, 2] + 0.5 * x[, 3]+1*(1.5+1.5*x[,2]-1.5*x[,3])))
  S0 <- pbinom(qceiling, size=1, prob = expit(-0.5+0.5 * x[, 2] + 0.5 * x[, 3]+0*(1.5+1.5*x[,2]-1.5*x[,3])), lower.tail = FALSE) + lambda * dbinom(qceiling, size=1, prob = expit(-0.5+0.5 * x[, 2] + 0.5 * x[, 3]+0*(1.5+1.5*x[,2]-1.5*x[,3])))

  return(as.numeric((S1 - S0) > 0))
}