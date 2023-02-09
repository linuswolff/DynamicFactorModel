#' A function to estimate a Dynamic Factor Model and use it to perform pseudo out of sample forecasts for the desired variable and the desired periods.
#'
#' @param X A matrix or data frame
#' @param variable The variable to be forecast
#' @param r The number of factors to be estimated
#' @param p The lag length of the ADL model
#' @param fcst_for_periods The periods to be forecast
#' @return The one step ahead forecasts and the real values of the supplied variable
#' @export 
loop_over_forecastDFM <- function(X, variable, r, p, fcst_for_periods) {
  
  fcsts_and_real <- cbind(numeric(length(fcst_for_periods)), X[fcst_for_periods, variable])
   
  for(i in 1:length(fcst_for_periods)) {
    fcsts_and_real[i, 1] <- forecastDFM(X = X[1:(fcst_for_periods[i] - 1),], variable = variable, r = r, p = p)
  } 
  return(list(fcsts_and_real = fcsts_and_real))
}
