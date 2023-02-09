#' A function to estimate a Dynamic Factor Model and perform forecasts with the estimated model
#'
#' @param X A matrix or data frame
#' @param variable The variable to be forecast
#' @param r The number of factors to be estimated
#' @param p The lag length of the ADL model
#' @return The one step ahead forecast of the supplied variable
#' @export 

forecastDFM <- function(X, variable, r, p) {
  library(dplyr)
  library(vars)
  
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
  
  factorVAR <- vars::VAR(F_hat, ic = "AIC", type = "none")
  fcst_F_1step <- predict(factorVAR, n.ahead = 1)
  fcst_Fplus1 <- c()
  for (i in 1:r) {
    fcst_Fplus1[i] <- fcst_F_1step$fcst[[i]][1]
  }
  
  # 3. combine into forecasting equation
  
  X_i_Tplus1 <- sum(coef(ADL.fit) * c(1, XwithLags[dim(XwithLags)[1], -1], fcst_Fplus1))
  
  X_i_Tplus1
  
}