#' A function to estimate the optimal number of factors for a Dynamic Factor Model.
#'
#' @param x A matrix or data frame
#' @param rmax The maximum number of factors to be considered
#' @param type Which of the three criterion from Bai and Ng (2002) should be used. (default = 2)
#' @return The optimal number of factors
#' @export 

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