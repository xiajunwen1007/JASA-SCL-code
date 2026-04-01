### define an exponential kernel
compute_kernel_matrix <- function(X, Y = NULL, sigma = 0.5, gamma = 0.1) {
  X <- as.matrix(X)
  if (is.null(Y)) Y <- X else Y <- as.matrix(Y)
  
  xx <- rowSums(X^2)
  yy <- rowSums(Y^2)
  sq_dist <- pmax(-2 * tcrossprod(X, Y) + outer(xx, yy, "+"), 0)
  
  exp(-gamma * sq_dist^sigma)
}

### Successive classification learning method for different exponential kernels when outcomes are count data.
#' @param data A list containing the outcome `y`, treatment indicator `a`, and covariates `x`.
#' @param tau The quantile level for which the quantile optimal treatment regime is to be estimated.
#' @param kappa A small positive number to control the precision of the binary search procedure.
#' @param epsilon A small positive number to control the convergence criterion of the binary search procedure.
#' @param mean.a The estimated propensity score pi(A=1|X).
#' @param f.y.ax1 The conditional survival model for treatment group.
#' @param f.y.ax0 The conditional survival model for control group. By combining with the f.y.ax1, we can construct the conditional negative binomial model to predict conditional survival function.
#' @param kernel The kernel to construct quantile optimal treatment regimes.
#' @param min The minimum value for the binary search procedure.
#' @param max The maximum value for the binary search procedure.
#' @param h The bandwidth for the soft regime Phi
#' @return A list containing opt.reg and q, representing the estimated quantile OTR and estimated quantile, respectively.
SCL.count.exponential <- function(data, tau, kappa=0, epsilon=0, mean.a, f.y.ax1, f.y.ax0,  kernel=0.5, min=0, max=0, h=0) {
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

  costs <- c(0.025/2, 0.05/2, 0.1/2, 1/2)
  gammas <- c(0.025/2, 0.05/2, 0.1/2, 1/2)
  K_cache <- lapply(gammas, function(gamma) {
    compute_kernel_matrix(x[,-1], sigma = kernel, gamma = gamma)
  })

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

    ### prepare data
    newa <- sign(psi1 - psi0)
    w    <- abs(psi1 - psi0)
    data.regime <- data.frame(y = factor(newa), x[, -1])

    ### fix folds once outside all loops
    folds <- sample(rep(1:10, length.out = nrow(data.regime)))

    ### parameter grid
    cv_errors_cg <- matrix(NA, nrow = length(costs), ncol = length(gammas),
                            dimnames = list(paste0("cost=", costs),
                                            paste0("gamma=", gammas)))

    ### CV loop — kernel precomputed once per gamma
    for (gi in seq_along(gammas)) {
      gamma <- gammas[gi]
      
      # compute full n x n kernel matrix ONCE per gamma
      K_full <- K_cache[[gi]]
      
      for (ci in seq_along(costs)) {
        cost <- costs[ci]
        cv_errors <- numeric(10)
        
        for (k in 1:10) {
          tr <- folds != k
          va <- folds == k
          
          # subset kernel matrix by fold — no recomputation needed
          K_train <- K_full[tr, tr, drop = FALSE]
          K_test  <- K_full[va, tr, drop = FALSE]
          
          opt.reg <- wsvm(K_train,
                          data.regime$y[tr],
                          cost   = cost,
                          weight = w[tr],
                          kernel = "precomputed")
          
          pred <- predict(opt.reg, K_test)
          # regime <- as.numeric(as.character(pred))
          cv_errors[k] <- sum(w[va] * (pred != data.regime$y[va]))
        }
        
        cv_errors_cg[ci, gi] <- mean(cv_errors)
      }
    }

    ### select best parameters
    best_params <- which(cv_errors_cg == min(cv_errors_cg), arr.ind = TRUE)
    best_cost   <- costs[best_params[1, 1]]
    best_gamma  <- gammas[best_params[1, 2]]

    ### fit final model on full data with best parameters
    K_final <- K_cache[[which(gammas == best_gamma)]]

    opt.reg <- wsvm(K_final,
                    data.regime$y,
                    cost   = best_cost,
                    weight = w,
                    kernel = "precomputed")

    ### calculate the optimal treatment regime
    dv  <- attr(predict(opt.reg, K_final, decision.values = TRUE), "decision.values")
    regime  <- Phi(newa[1] * dv, h)

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

  return(list(opt.reg=opt.reg, q=mid, best_gamma=best_gamma, X_train = x[,-1], kernel=kernel))
}

### Successive classification learning method for different exponential kernels with beta distribution for smoothing when outcomes are count data.
### Shape1 and shape2 are the parameters for the beta distribution. Other definitions of input are identical to that of SCL.count.exponential.
### The output is identical to SCL.count.exponential.
SCL.count.exponential.beta <- function(data, tau, kappa=0, epsilon=0, mean.a, f.y.ax1, f.y.ax0, shape1, shape2,  kernel=0.5, min=0, max=0, h=0) {
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

  costs <- c(0.025/2, 0.05/2, 0.1/2, 1/2)
  gammas <- c(0.025/2, 0.05/2, 0.1/2, 1/2)
  K_cache <- lapply(gammas, function(gamma) {
    compute_kernel_matrix(x[,-1], sigma = kernel, gamma = gamma)
  })

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

    ### prepare data
    newa <- sign(psi1 - psi0)
    w    <- abs(psi1 - psi0)
    data.regime <- data.frame(y = factor(newa), x[, -1])

    ### fix folds once outside all loops
    folds <- sample(rep(1:10, length.out = nrow(data.regime)))

    ### parameter grid
    cv_errors_cg <- matrix(NA, nrow = length(costs), ncol = length(gammas),
                            dimnames = list(paste0("cost=", costs),
                                            paste0("gamma=", gammas)))

    ### CV loop — kernel precomputed once per gamma
    for (gi in seq_along(gammas)) {
      gamma <- gammas[gi]
      
      # compute full n x n kernel matrix ONCE per gamma
      K_full <- K_cache[[gi]]
      
      for (ci in seq_along(costs)) {
        cost <- costs[ci]
        cv_errors <- numeric(10)
        
        for (k in 1:10) {
          tr <- folds != k
          va <- folds == k
          
          # subset kernel matrix by fold — no recomputation needed
          K_train <- K_full[tr, tr, drop = FALSE]
          K_test  <- K_full[va, tr, drop = FALSE]
          
          opt.reg <- wsvm(K_train,
                          data.regime$y[tr],
                          cost   = cost,
                          weight = w[tr],
                          kernel = "precomputed")
          
          pred <- predict(opt.reg, K_test)
          # regime <- as.numeric(as.character(pred))
          cv_errors[k] <- sum(w[va] * (pred != data.regime$y[va]))
        }
        
        cv_errors_cg[ci, gi] <- mean(cv_errors)
      }
    }

    ### select best parameters
    best_params <- which(cv_errors_cg == min(cv_errors_cg), arr.ind = TRUE)
    best_cost   <- costs[best_params[1, 1]]
    best_gamma  <- gammas[best_params[1, 2]]

    ### fit final model on full data with best parameters
    K_final <- K_cache[[which(gammas == best_gamma)]]

    opt.reg <- wsvm(K_final,
                    data.regime$y,
                    cost   = best_cost,
                    weight = w,
                    kernel = "precomputed")

    ### calculate the optimal treatment regime
    dv  <- attr(predict(opt.reg, K_final, decision.values = TRUE), "decision.values")
    regime  <- Phi(newa[1] * dv, h)

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

  return(list(opt.reg=opt.reg, q=mid, best_gamma=best_gamma, X_train = x[,-1], kernel=kernel))
}

expit <- function(x) {
  return(exp(x) / (1 + exp(x)))
}

### A function to calculate the difference between the optimal survival function at q and 1-tau. 
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

### A function to calculate the optimal treatment that maximizes the survival function at point q.
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

### A function to calculate the difference between the optimal survival function at q and 1-tau. 
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

### A function to calculate the optimal treatment that maximizes the survival function at point q.
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

### A prediction function for SCL objects
### Input the SCL object and new data to predict.
### Output the predicted values.
predict.SCL <- function(model, x.test, chunk_size = 10000) {
  X_test  <- as.matrix(x.test)
  X_train <- as.matrix(model$X_train)
  n_test  <- nrow(X_test)
  
  preds <- integer(n_test)
  
  for (start in seq(1, n_test, by = chunk_size)) {
    end     <- min(start + chunk_size - 1, n_test)
    X_chunk <- X_test[start:end, , drop = FALSE]
    
    K_chunk  <- compute_kernel_matrix(X_chunk, X_train, sigma = model$kernel, gamma = model$best_gamma)
    
    preds[start:end] <- as.integer(predict(model$opt.reg, K_chunk))
    
    rm(X_chunk, K_chunk); gc()
  }
  
  preds
}