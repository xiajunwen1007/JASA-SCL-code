### function to perform grid search method for estimating quantile OTRs
### data: a list containing x, a, y; tau: quantile level.
grid_serach <- function(data, tau){
  # extract x, a, y
  x <- data$x
  a <- data$a
  y <- data$y

  # We separately calculate the quantile for beta0 = 1, beta0 = -1, and beta0 = 0 for beta=df
  beta <- expand.grid(beta1 = seq(0, 6, length.out = 15),
                  beta2 = seq(-1, 1, length.out = 5),
                  beta3=seq(-1,1,length.out=5),
                  beta4=seq(-1,1,length.out=5),
                  beta5=seq(-1,1,length.out=5))

  quantile.ipw.1 <- function(i){
   quant_est(t(cbind(1,beta[i, ])), x, y, a, rep(0.5,n), tau = tau )
  }
  quant1 <- sapply(1:nrow(beta),quantile.ipw.1)

  quantile.ipw.minus1 <- function(i){
   quant_est(t(cbind(-1,beta[i, ])), x, y, a, rep(0.5,n), tau=tau)
  }
  quantminus1 <- sapply(1:nrow(beta),quantile.ipw.minus1)

  quantile.ipw.0 <- function(i){
   quant_est(t(cbind(0,beta[i, ])), x, y, a, rep(0.5,n), tau=tau)
  }
  quant0 <- sapply(1:nrow(beta),quantile.ipw.0)

  # calculate the maximum quantile among these regimes
  maxi <- max(c(quant1,quantminus1,quant0))

  # calculate the coefficients 
  if(max(quant1)>max(c(quantminus1,quant0))){
    coef <- t(cbind(1,beta[which(c(quant1)==max(quant1)),][1,]))
  }else if (max(quantminus1)>max(quant0)){
    coef <- t(cbind(-1,beta[which(c(quantminus1)==max(quantminus1)),][1,]))
  }else{
    coef <- t(cbind(0,beta[which(c(quant0)==max(quant0)),][1,]))
  }

  return(list(coef=coef, hatQ=maxi))
 }

### Function to calculate the difference between the optimal survival function at q and 1-tau. 
opt.q <- function(q, tau) {
  opt.a <- opt.reg(q)
  return(mean(pnorm(q, mean = eval(parse(text = base.test)) + opt.a * (x.test[, 2] + x.test[, 3]+x.test[,4]+x.test[,5]+x.test[,6]), sd = 2 * opt.a + 1, lower.tail = FALSE)) - 1 + tau)
}

### Function to calculate the optimal treatment that maximizes the survival function at point q.
opt.reg <- function(q) {
  return(as.numeric((q - eval(parse(text = base.test))) / 1 - ((q - eval(parse(text = base.test)) - 1 * x.test[, 2] -1 * x.test[, 3]-1 * x.test[,4]-1 *x.test[,5]-1*x.test[,6]) / 3) > 0))
}