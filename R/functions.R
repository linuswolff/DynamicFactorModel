# Hello, world!
#
# This is an example function named 'hello' 
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

forecastDFM <- function(X, variable, r, p, max.lag.var) {
  library(dplyr)
  X_mat <- as.matrix(X)
  n <- ncol(X_mat)
  t <- nrow(X_mat)
  
  ## DFM estimation
  
  # n x n sample cov matrix
  cov_mat_hat <- (t(X_mat) %*% X_mat) / t
  
  # n x r matrix of normalized eigenvectors
  V_hat <- eigen(cov_mat_hat)$vectors[, 1:r]
  
  # estimates of Lambda, Factors and Common Component
  L_hat <- V_hat * sqrt(n)
  F_hat <- (X_mat %*% V_hat) / sqrt(n)
  C_hat <- F_hat %*% t(L_hat)
  idio_hat <- X_mat - C_hat
  
  ## forecasting
  # 1. AR model
  
  X_i <- X_mat[, variable]
  len <- length(X_i)
  XwithLags <- matrix(c(X_i, numeric(len * p)), len, p + 1, byrow = FALSE)
  for (i in 1:p) {
    XwithLags[, i + 1] <- dplyr::lag(X_i, i)
  }
  XwithLags <- XwithLags[(p + 1):len, ]
  
  
  ADL.fit <- lm(XwithLags[, 1] ~ XwithLags[, -1] + F_hat[(p + 1):len, 1:r])
  
  # 2. VAR model for factors
  
  factorVAR <- VAR(F_hat, ic = "AIC", type = "none")
  fcst_F_1step <- predict(factorVAR, n.ahead = 1)
  fcst_Fplus1 <- c()
  for (i in 1:r) {
    fcst_Fplus1[i] <- fcst_F_1step$fcst[[i]][1]
  }
  
  # 3. combine into forecasting equation
  
  X_i_Tplus1 <- sum(coef(ADL.fit) * c(1, XwithLags[dim(XwithLags)[1], -1], fcst_Fplus1))
  
  X_i_Tplus1
  
}

loop_over_forecastDFM <- function(X, variable, r, p, max.lag.var, fcst_for_periods) {
  
  fcsts_and_real <- cbind(numeric(length(fcst_for_periods)), X[fcst_for_periods, variable])
  
  for(i in 1:length(fcst_for_periods)) {
    fcsts_and_real[i, 1] <- forecastDFM(X = X[1:(fcst_for_periods[i] - 1),], variable = variable, r = r, p = p, max.lag.var = max.lag.var)
  } 
  return(list(fcsts_and_real = fcsts_and_real))
}

ICfactors_mine <- function (x, rmax, type = 2) {
  x <- as.matrix(x)
  Mx <- colMeans(x)
  Wx <- apply(x, MARGIN = 2, FUN = sd)
  
  for (i in 1:ncol(x)) {
    x[, i] <- (x[, i] - Mx[i])/Wx[i]
  }
  
  TT <- nrow(x)
  N <- ncol(x)
  eigen <- eigen(cov(x))
  result <- numeric(rmax)
  
  if (type == 1) {
    for (r in 1:rmax) {
      v <- eigen$vectors[, 1:r]
      factors <- x %*% v
      V <- sum(diag(t(x - factors %*% t(v)) %*% (x - factors %*% t(v)))/(nrow(x) * ncol(x)))
      result[r] <- log(V) + r * (N + TT)/(N * TT) * log(N * TT/(N + TT))
    }
    
    #plot(result, main = "ICR1", xlab = "Number of fators", ylab = "Index")
    #points(which.min(result), result[which.min(result)], pch = 19, col = "red")
    
    list(r_star = which.min(result)[1], IC = result)
  }
  else if (type == 2) {
    for (r in 1:rmax) {
      v <- eigen$vectors[, 1:r]
      factors <- x %*% v
      V <- sum(diag(t(x - factors %*% t(v)) %*% (x - factors %*% t(v)))/(nrow(x) * ncol(x)))
      result[r] <- log(V) + r * (N + TT)/(N * TT) * log(min(N, TT))
    }
    
    #plot(result, main = "ICR2", xlab = "Number of fators", ylab = "Index")
    #points(which.min(result), result[which.min(result)], pch = 19, col = "red")
    
    list(r_star = which.min(result)[1], IC = result)
  }
  else if (type == 3) {
    for (r in 1:rmax) {
      v <- eigen$vectors[, 1:r]
      factors <- x %*% v
      V <- sum(diag(t(x - factors %*% t(v)) %*% (x - factors %*% t(v)))/(nrow(x) * ncol(x)))
      result[r] <- log(V) + r * (log(min(N, TT))/min(N, TT))
    }
    
    #plot(result, main = "ICR3", xlab = "Number of fators", ylab = "Index")
    #points(which.min(result), result[which.min(result)], pch = 19, col = "red")
    
    list(r_star = which.min(result)[1], IC = result)
  }
}