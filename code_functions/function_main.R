### Successive classification learning method for estimating quantile optimal treatment regimes when outcomes are continuous.
#' @param data A list containing the outcome `y`, treatment indicator `a`, and covariates `x`.
#' @param tau The quantile level for which the quantile optimal treatment regime is to be estimated.
#' @param kappa A small positive number to control the precision of the binary search procedure.
#' @param epsilon A small positive number to control the convergence criterion of the binary search procedure.
#' @param mean.a The estimated propensity score pi(A=1|X).
#' @param f.y.ax The conditional mean model f.y.ax.
#' @param f.r.ax The conditional variance model f.r.ax. By combining with the f.y.ax, we can construct the conditional gaussian model to predict conditional survival function.
#' @param kernel The kernel type to construct quantile optimal treatment regimes, which can be "gaussian" and "linear". Default is "linear".
#' @param min The minimum value for the binary search procedure.
#' @param max The maximum value for the binary search procedure.
#' @param h The bandwidth for the soft regime Phi
#' @return A list containing opt.reg and q, representing the estimated quantile OTR and estimated quantile, respectively.
SCL.continuous <- function(data, tau, kappa=0, epsilon=0, mean.a, f.y.ax, f.r.ax,  kernel="linear", min=0, max=0, h=0) {
  # extract data components
  y <- data$y
  a <- data$a
  x <- data$x
  n <- nrow(x)

  # initialize the binary search procedure
  if (min == 0 && max == 0) {
    min <- min(y)
    max <- max(y)
  }
  mid <- (max + min) / 2

  # parameters will be used in the binary search procedure
  # the smoothed regime function
  Phi <-function(x,h) pnorm(x/h)

  # the bandwidth for the smoothed regime function
  if (h==0){
     h <- 0.2/log(n)
  }
 
  # kappa
  if (kappa == 0) {
    kappa <- sd(y) / sqrt(n) / 6
  }
  # epsilon
  if (epsilon == 0) {
    epsilon <- 1 / sqrt(n) / 2
  }

  # the binary search procedure
  for (iter in 1:100) {
    newx1 <- data.frame(y=y,x=x[,-1], a=1)
    newx0 <- data.frame(y=y,x=x[,-1], a=0)

    # conditional survival function at mid
    S1 <- pnorm(mid, mean = predict(f.y.ax, newx1), sd = sqrt(predict(f.r.ax, newx1, type = "response")), lower.tail = FALSE)
    S0 <- pnorm(mid, mean = predict(f.y.ax, newx0), sd = sqrt(predict(f.r.ax, newx0, type = "response")), lower.tail = FALSE)
    
    # psi function values to achieve the doubly robust property
    psi1 <- S1 + (as.numeric(y > mid) - S1) / mean.a * a
    psi0 <- S0 + (as.numeric(y > mid) - S0) / (1 - mean.a) * (1 - a)

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
      tn <- tune_wsvm(y ~ ., data = data.regime, kernel = "linear", ranges = list(cost = c(0.025 / 2, 0.05 / 2, 0.1 / 2, 1 / 2)), weight = w)
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

### Successive classification learning method for estimating quantile optimal treatment regimes when outcomes are count data.
#' @param data A list containing the outcome `y`, treatment indicator `a`, and covariates `x`.
#' @param tau The quantile level for which the quantile optimal treatment regime is to be estimated.
#' @param kappa A small positive number to control the precision of the binary search procedure.
#' @param epsilon A small positive number to control the convergence criterion of the binary search procedure.
#' @param mean.a The estimated propensity score pi(A=1|X).
#' @param f.y.ax1 The conditional survival model for treatment group.
#' @param f.y.ax0 The conditional survival model for control group. By combining with the f.y.ax1, we can construct the conditional negative binomial model to predict conditional survival function.
#' @param kernel The kernel type to construct quantile optimal treatment regimes, which can be "gaussian" and "linear". Default is "linear".
#' @param min The minimum value for the binary search procedure.
#' @param max The maximum value for the binary search procedure.
#' @param h The bandwidth for the soft regime Phi
#' @return A list containing opt.reg and q, representing the estimated quantile OTR and estimated quantile, respectively.
SCL.count <- function(data, tau, kappa=0, epsilon=0, mean.a, f.y.ax1, f.y.ax0,  kernel="linear", min=0, max=0, h=0) {
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
    # calculate the conditional smoothed survival function at mid
    qfloor <- floor(mid)
    qceiling <- ceiling(mid)
    if (qfloor == qceiling) {
      lambda <- 0
    } else {
      lambda <- (mid - qceiling) / (qfloor - qceiling)
    }

    # conditional smoothed survival function at mid
    S1 <- pnbinom(qceiling, mu = predict(f.y.ax1, data.frame(x=x[,-1]), type = "response"), size = f.y.ax1$theta, lower.tail = FALSE) + lambda * dnbinom(qceiling, mu = predict(f.y.ax1, data.frame(x=x[,-1]), type = "response"), size = f.y.ax1$theta)

    S0 <- pnbinom(qceiling, mu = predict(f.y.ax0, data.frame(x=x[,-1]), type = "response"), size = f.y.ax0$theta, lower.tail = FALSE) + lambda * dnbinom(qceiling, mu = predict(f.y.ax0, data.frame(x=x[,-1]), type = "response"), size = f.y.ax0$theta)
    
    # smoothed psi function values to achieve the doubly robust property
    psi1 <- S1 + (as.numeric(y > qceiling) + lambda * (as.numeric(y == qceiling)) - S1) / mean.a * a
    psi0 <- S0 + (as.numeric(y > qceiling) + lambda * (as.numeric(y == qceiling)) - S0) / (1 - mean.a) * (1 - a)

    # calculate the new treatment indicator and the weights
    newa <- sign(psi1 - psi0)
    w <- abs(psi1 - psi0)
    data.regime <- data.frame(y = factor(newa), x = x[,-1])

    # if the kernel is "gaussian", use the gaussian kernel, otherwise use the linear kernel
    if (kernel == "gaussian") {
      # tune the cost parameter for the weighted SVM
       tn <- tune_wsvm(y ~ ., data = data.regime, kernel = "radial", ranges = list(cost = c(0.025 / 2, 0.05 / 2, 0.1 / 2, 1 / 2), gamma = c(0.025 / 2, 0.05 / 2, 0.1 / 2, 1 / 2)), weight = w)
    
      # fit the weighted SVM model
      opt.reg <- wsvm(y ~ ., data = data.regime, kernel = "radial", cost = tn$best.parameters$cost, gamma = tn$best.parameters$gamma, weight = w)

      # calculate the optimal treatment regime
      regime = Phi(newa[1]*attr(predict(opt.reg, data.regime[,-1], decision.values=TRUE),"decision.values"),h)
      
    } else {
      # use the linear kernel
      tn <- tune_wsvm(y ~ ., data = data.regime, kernel = "linear", ranges = list(cost = c(0.025 / 2, 0.05 / 2, 0.1 / 2, 1 / 2)), weight = w)
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

### Successive classification learning method with the hard regime for estimating quantile optimal treatment regimes when outcomes are continuous.
### The definition of input is identical to that of SCL.continuous.
### The output is identical to SCL.continuous.
SCL.continuous.hard <- function(data, tau, kappa=0, epsilon=0, mean.a, f.y.ax, f.r.ax,  kernel="linear", min=0, max=0, h=0) {
  # extract data components
  y <- data$y
  a <- data$a
  x <- data$x
  n <- nrow(x)

  # initialize the binary search procedure
  if (min == 0 && max == 0) {
    min <- min(y)
    max <- max(y)
  }
  mid <- (max + min) / 2

  # parameters will be used in the binary search procedure 
  # kappa
  if (kappa == 0) {
    kappa <- sd(y) / sqrt(n) / 6
  }
  # epsilon
  if (epsilon == 0) {
    epsilon <- 1 / sqrt(n) / 2
  }

  # the binary search procedure
  for (iter in 1:100) {
    newx1 <- data.frame(y=y, x=x[,-1], a=1)
    newx0 <- data.frame(y=y, x=x[,-1], a=0)

    # conditional survival function at mid
    S1 <- pnorm(mid, mean = predict(f.y.ax, newx1), sd = sqrt(predict(f.r.ax, newx1, type = "response")), lower.tail = FALSE)
    S0 <- pnorm(mid, mean = predict(f.y.ax, newx0), sd = sqrt(predict(f.r.ax, newx0, type = "response")), lower.tail = FALSE)
    
    # psi function values to achieve the doubly robust property
    psi1 <- S1 + (as.numeric(y > mid) - S1) / mean.a * a
    psi0 <- S0 + (as.numeric(y > mid) - S0) / (1 - mean.a) * (1 - a)

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
      regime = as.numeric(predict(opt.reg, data.regime[,-1]))-1
    } else {
      # use the linear kernel
      tn <- tune_wsvm(y ~ ., data = data.regime, kernel = "linear", ranges = list(cost = c(0.025 / 2, 0.05 / 2, 0.1 / 2, 1 / 2)), weight = w)
      opt.reg <- wsvm(y ~ ., data = data.regime, kernel = "linear", cost = tn$best.parameters$cost, weight = w)
      regime = as.numeric(predict(opt.reg, data.regime[,-1]))-1
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

### Successive classification learning method with the hard regime for estimating quantile optimal treatment regimes when outcomes are count data.
### The definition of input is identical to that of SCL.count.
### The output is identical to SCL.count.
SCL.count.hard <- function(data, tau, kappa=0, epsilon=0, mean.a, f.y.ax1, f.y.ax0,  kernel="linear", min=0, max=0, h=0) {
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
    # calculate the conditional smoothed survival function at mid
    qfloor <- floor(mid)
    qceiling <- ceiling(mid)
    if (qfloor == qceiling) {
      lambda <- 0
    } else {
      lambda <- (mid - qceiling) / (qfloor - qceiling)
    }

    # conditional smoothed survival function at mid
    S1 <- pnbinom(qceiling, mu = predict(f.y.ax1, data.frame(x=x[,-1]), type = "response"), size = f.y.ax1$theta, lower.tail = FALSE) + lambda * dnbinom(qceiling, mu = predict(f.y.ax1, data.frame(x=x[,-1]), type = "response"), size = f.y.ax1$theta)

    S0 <- pnbinom(qceiling, mu = predict(f.y.ax0, data.frame(x=x[,-1]), type = "response"), size = f.y.ax0$theta, lower.tail = FALSE) + lambda * dnbinom(qceiling, mu = predict(f.y.ax0, data.frame(x=x[,-1]), type = "response"), size = f.y.ax0$theta)
    
    # smoothed psi function values to achieve the doubly robust property
    psi1 <- S1 + (as.numeric(y > qceiling) + lambda * (as.numeric(y == qceiling)) - S1) / mean.a * a
    psi0 <- S0 + (as.numeric(y > qceiling) + lambda * (as.numeric(y == qceiling)) - S0) / (1 - mean.a) * (1 - a)

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
      regime = as.numeric(predict(opt.reg, data.regime[,-1]))-1
    } else {
      # use the linear kernel
      tn <- tune_wsvm(y ~ ., data = data.regime, kernel = "linear", ranges = list(cost = c(0.025 / 2, 0.05 / 2, 0.1 / 2, 1 / 2)), weight = w)
      opt.reg <- wsvm(y ~ ., data = data.regime, kernel = "linear", cost = tn$best.parameters$cost, weight = w)
      regime = as.numeric(predict(opt.reg, data.regime[,-1]))-1
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

### Successive classification learning method with beta distribution for smoothing when outcomes are count data.
### Shape1 and shape2 are the parameters for the beta distribution. Other definitions of input are identical to that of SCL.count.
### The output is identical to SCL.count.
SCL.count.beta <- function(data, tau, kappa=0, epsilon=0, mean.a, f.y.ax1, f.y.ax0,  kernel="linear", shape1, shape2, min=0, max=0, h=0) {
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
    # calculate the conditional smoothed survival function at mid
    qfloor <- floor(mid)
    qceiling <- ceiling(mid)
    if (qfloor == qceiling) {
      lambda <- 0
    } else {
      lambda <- (mid - qceiling) / (qfloor - qceiling)
      lambda <- pbeta(lambda, shape1, shape2)
    }

    # conditional smoothed survival function at mid
    S1 <- pnbinom(qceiling, mu = predict(f.y.ax1, data.frame(x=x[,-1]), type = "response"), size = f.y.ax1$theta, lower.tail = FALSE) + lambda * dnbinom(qceiling, mu = predict(f.y.ax1, data.frame(x=x[,-1]), type = "response"), size = f.y.ax1$theta)

    S0 <- pnbinom(qceiling, mu = predict(f.y.ax0, data.frame(x=x[,-1]), type = "response"), size = f.y.ax0$theta, lower.tail = FALSE) + lambda * dnbinom(qceiling, mu = predict(f.y.ax0, data.frame(x=x[,-1]), type = "response"), size = f.y.ax0$theta)
    
    # smoothed psi function values to achieve the doubly robust property
    psi1 <- S1 + (as.numeric(y > qceiling) + lambda * (as.numeric(y == qceiling)) - S1) / mean.a * a
    psi0 <- S0 + (as.numeric(y > qceiling) + lambda * (as.numeric(y == qceiling)) - S0) / (1 - mean.a) * (1 - a)

    # calculate the new treatment indicator and the weights
    newa <- sign(psi1 - psi0)
    w <- abs(psi1 - psi0)
    data.regime <- data.frame(y = factor(newa), x = x[,-1])

    # if the kernel is "gaussian", use the gaussian kernel, otherwise use the linear kernel
    if (kernel == "gaussian") {
      # tune the cost parameter for the weighted SVM
       tn <- tune_wsvm(y ~ ., data = data.regime, kernel = "radial", ranges = list(cost = c(0.025 / 2, 0.05 / 2, 0.1 / 2, 1 / 2), gamma = c(0.025 / 2, 0.05 / 2, 0.1 / 2, 1 / 2)), weight = w)
    
      # fit the weighted SVM model
      opt.reg <- wsvm(y ~ ., data = data.regime, kernel = "radial", cost = tn$best.parameters$cost, gamma = tn$best.parameters$gamma, weight = w)

      # calculate the optimal treatment regime
      regime = Phi(newa[1]*attr(predict(opt.reg, data.regime[,-1], decision.values=TRUE),"decision.values"),h)
      
    } else {
      # use the linear kernel
      tn <- tune_wsvm(y ~ ., data = data.regime, kernel = "linear", ranges = list(cost = c(0.025 / 2, 0.05 / 2, 0.1 / 2, 1 / 2)), weight = w)
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

### QIQ-learning for estimating quantile optimal treatment regimes
### The definition of the input is identical to that of SCL.continuous. 
### The output is identical to SCL.continuous.
QIQ.learning.continuous <- function(data, tau, f.y.ax, f.r.ax) {
  # extract data components
  y <- data$y
  a <- data$a
  x <- data$x

  # construct a function opt such that for any q we can obtain the difference between the estimated optimal survival function at q and 1-tau
  newx1 <- data.frame(y=y, x=x[,-1], a=1)
  newx0 <- data.frame(y=y, x=x[,-1], a=0)
  opt <- function(q) {
    # estimated optimal treatments that maximizes survival function at q
    opt.a <- as.numeric((q - predict(f.y.ax, newx0)) / sqrt(predict(f.r.ax, newx0, type = "response")) - (q - predict(f.y.ax, newx1)) / sqrt(predict(f.r.ax, newx1, type = "response")) > 0)

    opt.data <- data.frame(y=y, x=x[,-1], a=opt.a)

    return(mean(pnorm(q, mean = predict(f.y.ax, opt.data), sd = sqrt(predict(f.r.ax, opt.data, type = "response")), lower.tail = FALSE)) - 1 + tau)
  }

  # estimate the optimal quantile
  opt.q <- uniroot(opt, interval = c(min(y), max(y)))$root

  # construct an object of opt.reg
  opt.reg <- list(
    f.y.ax = f.y.ax,
    f.r.ax = f.r.ax,
    opt.q = opt.q
  )

  # give it an S3 class
  class(opt.reg) <- "QIQlearningcontinuous"

  return(list(opt.reg=opt.reg, q = opt.q))
}

### QIQ-learning for estimating quantile optimal treatment regimes
### The definition of the input is identical to of SCL.count. 
### The output is identical to SCL.count
QIQ.learning.count <- function(data, tau, f.y.ax1, f.y.ax0) {
  # extract data components
  y <- data$y
  a <- data$a
  x <- data$x

  # construct a function opt such that for any q we can obtain the difference between the estimated optimal survival function at q and 1-tau
  opt <- function(q) {
    qfloor <- floor(q)
    qceiling <- ceiling(q)

    if (qfloor == qceiling) {
      lambda <- 0
    } else {
      lambda <- (q - qceiling) / (qfloor - qceiling)
    }

    # estimated optimal treatments that maximizes survival function at q
    opt.a <- as.numeric(pnbinom(qceiling, mu = predict(f.y.ax0, data.frame(x=x[, -1]), type = "response"), size = f.y.ax0$theta, lower.tail = FALSE)+lambda*dnbinom(qceiling, mu=predict(f.y.ax0, data.frame(x=x[, -1]), type = "response"), size = f.y.ax0$theta) - pnbinom(qceiling, mu = predict(f.y.ax1, data.frame(x=x[, -1]), type = "response"), size = f.y.ax1$theta, lower.tail = FALSE)- lambda*dnbinom(qceiling, mu = predict(f.y.ax1, data.frame(x=x[, -1]), type = "response"), size = f.y.ax1$theta)< 0)
    
    return(mean(pnbinom(qceiling, mu = opt.a * predict(f.y.ax1, data.frame(x=x[, -1]), type = "response") + (1 - opt.a) * predict(f.y.ax0, data.frame(x=x[, -1]), type = "response"), size = f.y.ax1$theta * opt.a + f.y.ax0$theta * (1 - opt.a), lower.tail = FALSE)+lambda*dnbinom(qceiling, mu = opt.a * predict(f.y.ax1, data.frame(x=x[, -1]), type = "response") + (1 - opt.a) * predict(f.y.ax0, data.frame(x=x[, -1]), type = "response"), size = f.y.ax1$theta * opt.a + f.y.ax0$theta * (1 - opt.a))) - 1 + tau)
  }

  # estimate the optimal quantile
  opt.q <- uniroot(opt, c(-1, max(y)))$root
  

  # construct an object of opt.reg
  opt.reg <- list(
    f.y.ax1 = f.y.ax1,
    f.y.ax0 = f.y.ax0,
    opt.q = opt.q
  )

  # give it an S3 class
  class(opt.reg) <- "QIQlearningcount"

  return(list(opt.reg=opt.reg, q = opt.q))
}

### Some utility functions
### A prediction function for the object opt.reg outputed from QIQ.learning.continuous function.
### Input the opt.reg outputed from QIQ.learning.continuous function and the covariates x.
### Output the predicted optimal treatment decisions.
predict.QIQlearningcontinuous <- function(object, x, ...) {
  newx1 <- data.frame(y=1, x=x, a=1)
  newx0 <- data.frame(y=1, x=x, a=0)
  as.numeric(
    (object$opt.q - predict(object$f.y.ax, newx0)) / sqrt(predict(object$f.r.ax, newx0, type = "response")) -
    (object$opt.q - predict(object$f.y.ax, newx1)) / sqrt(predict(object$f.r.ax, newx1, type = "response")) > 0
  )
}

### A prediction function for the object opt.reg outputed from QIQ.learning.count function.
### Input the opt.reg outputed from QIQ.learning.count function and the covariates x.
### Output the predicted optimal treatment decisions.
predict.QIQlearningcount <- function(object, x, ...) {
  opt.q <- object$opt.q
  qceiling <- ceiling(opt.q)
  qfloor <- floor(opt.q)
  if (qfloor == qceiling) {
      lambda <- 0
    } else {
      lambda <- (opt.q - qceiling) / (qfloor - qceiling)
    }

  as.numeric(pnbinom(qceiling, mu = predict(object$f.y.ax0, data.frame(x=x), type = "response"), size = object$f.y.ax0$theta, lower.tail = FALSE)+lambda*dnbinom(qceiling, mu=predict(object$f.y.ax0, data.frame(x=x), type = "response"), size = object$f.y.ax0$theta) - pnbinom(qceiling, mu = predict(object$f.y.ax1, data.frame(x=x), type = "response"), size = object$f.y.ax1$theta, lower.tail = FALSE)- lambda*dnbinom(qceiling, mu = predict(object$f.y.ax1, data.frame(x=x), type = "response"), size = object$f.y.ax1$theta)< 0)
}

### A prediction function for glmnetformula objects
### Input the glmnetformula object and new data to predict.
### Output the predicted values.
predict.glmnetformula <- function(object, new_data, type="response", ...) {
  X_new <- model.matrix(object$formula, new_data)[, -1]
  predict(object$fit, newx = X_new, type= type, ...)
}

### functions to calculate the optimal treatment that maximizes the survival function at point q.
opt.reg.case1 <- function(q) {
  return(as.numeric((q - eval(parse(text = case1.base)) ) / 1 - ((q - eval(parse(text = case1.base))- x[, 2] - x[, 3]) / 3) > 0))
}

opt.reg.case2 <- function(q){
 return(as.numeric((q - eval(parse(text = case2.base)))/exp(-2-0.5*x[,2]-0.5*x[,3])-((q - eval(parse(text = case2.base)) - x[,2]-x[,3])/exp(0.5*x[,2]+0.5*x[,3]))>0))
}

expit <- function(x) {
  return(exp(x) / (1 + exp(x)))
}

opt.reg.case3 <- function(q) {
  qfloor <- floor(q)
  qceiling <- ceiling(q)

  if (qfloor == qceiling) {
    lambda <- 0
  } else {
    lambda <- (q - qceiling) / (qfloor - qceiling)
  }

  S0 <- pnbinom(qceiling, prob = expit(eval(parse(text = case3.base))), lower.tail = FALSE, size = 1) + lambda * dnbinom(qceiling, prob = expit(eval(parse(text = case3.base))), size = 1)

  S1 <- pnbinom(qceiling, prob = expit(eval(parse(text = case3.base)) + 1 + x[, 2] + x[, 3]^2), lower.tail = FALSE, size = 3) + lambda * dnbinom(qceiling, prob = expit(eval(parse(text = case3.base)) + 1+ x[, 2] + x[, 3]^2), size = 3)

  return(as.numeric((S1 - S0) > 0))
}

# shape1 and shape2 are the parameters of the beta distribution used for smoothing
opt.reg.case3.beta <- function(q, shape1, shape2) {
  qfloor <- floor(q)
  qceiling <- ceiling(q)

  if (qfloor == qceiling) {
    lambda <- 0
  } else {
    lambda <- (q - qceiling) / (qfloor - qceiling)
    lambda <- pbeta(lambda, shape1, shape2)
  }

  S0 <- pnbinom(qceiling, prob = expit(eval(parse(text = case3.base))), lower.tail = FALSE, size = 1) + lambda * dnbinom(qceiling, prob = expit(eval(parse(text = case3.base))), size = 1)

  S1 <- pnbinom(qceiling, prob = expit(eval(parse(text = case3.base)) + 1 + x[, 2] + x[, 3]^2), lower.tail = FALSE, size = 3) + lambda * dnbinom(qceiling, prob = expit(eval(parse(text = case3.base)) + 1+ x[, 2] + x[, 3]^2), size = 3)

  return(as.numeric((S1 - S0) > 0))
}

### functions to calculate the difference between the optimal survival function at q and 1-tau. 
opt.q.case1 <- function(q, tau) {
  opt.a <- opt.reg.case1(q)
  return(mean(pnorm(q, mean = eval(parse(text = case1.base)) + opt.a * (x[, 2] + x[, 3]), sd = 2 * opt.a + 1, lower.tail = FALSE)) - 1 + tau)
}

opt.q.case2 <- function(q, tau){
 opt.a <- opt.reg.case2(q)
 return(mean(pnorm(q, mean=eval(parse(text = case2.base))+opt.a*(x[,2]+x[,3]),sd=exp(2*opt.a-2+0.5*(2*opt.a-1)*(x[,2]+x[,3])),lower.tail = FALSE))-1+tau)
}

opt.q.case3 <- function(q, tau) {
  opt.a <- opt.reg.case3(q)

  qfloor <- floor(q)
  qceiling <- ceiling(q)

  if (qfloor == qceiling) {
    lambda <- 0
  } else {
    lambda <- (q - qceiling) / (qfloor - qceiling)
  }
  return(mean(pnbinom(qceiling, prob = expit(eval(parse(text = case3.base)) + opt.a * (1 + x[, 2] + x[, 3]^2)), size = 2 * opt.a + 1, lower.tail = FALSE) + lambda * dnbinom(qceiling, prob = expit(eval(parse(text = case3.base)) + opt.a * (1 + x[, 2] + x[, 3]^2)), size = 2 * opt.a + 1)) - 1 + tau)
}

# Shape1 and shape2 are the parameters of the beta distribution used for smoothing
opt.q.case3.beta <- function(q, tau, shape1, shape2) {
  opt.a <- opt.reg.case3.beta(q, shape1, shape2)

  qfloor <- floor(q)
  qceiling <- ceiling(q)

  if (qfloor == qceiling) {
    lambda <- 0
  } else {
    lambda <- (q - qceiling) / (qfloor - qceiling)
    lambda <- pbeta(lambda, shape1, shape2)
  }

  return(mean(pnbinom(qceiling, prob = expit(eval(parse(text = case3.base)) + opt.a * (1 + x[, 2] + x[, 3]^2)), size = 2 * opt.a + 1, lower.tail = FALSE) + lambda * dnbinom(qceiling, prob = expit(eval(parse(text = case3.base)) + opt.a * (1 + x[, 2] + x[, 3]^2)), size = 2 * opt.a + 1)) - 1 + tau)
}

### A function to calculate estimated smooth quantile
hat.smooth.quantile <- function(y, a, mean.a, regime, tau) {
  sort.y <- sort(unique(y))
  surv <- sapply(sort.y, function(q) mean(as.numeric(y > q)*((a==regime)* (1/mean.a * a + 1/(1-mean.a)*(1-a)) )))

  qceiling.pos <- which.min(surv >= (1 - tau))

  if (qceiling.pos == "-Inf") {
    qceiling.pos <- 1
  }

  qfloor.pos <- qceiling.pos - 1

  lambda <- (1 - tau - surv[qceiling.pos]) / (ifelse(qfloor.pos == 0, 1, surv[qfloor.pos]) - surv[qceiling.pos])

  q <- lambda * sort.y[qfloor.pos] + (1 - lambda) * sort.y[qceiling.pos]

  return(q)
}

### A function to calculate true smooth quantile
smooth.quantile <- function(y, tau) {
  sort.y <- sort(unique(y))
  surv <- sapply(sort.y, function(q) mean(as.numeric(y > q)))

  qceiling.pos <- which.min(surv >= (1 - tau))

  if (qceiling.pos == "-Inf") {
    qceiling.pos <- 1
  }

  qfloor.pos <- qceiling.pos - 1

  lambda <- (1 - tau - surv[qceiling.pos]) / (ifelse(qfloor.pos == 0, 1, surv[qfloor.pos]) - surv[qceiling.pos])

  q <- lambda * sort.y[qfloor.pos] + (1 - lambda) * sort.y[qceiling.pos]

  return(q)
}