% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loop_over_forecastDFM.R
\name{loop_over_forecastDFM}
\alias{loop_over_forecastDFM}
\title{A function to estimate a Dynamic Factor Model and use it to perform pseudo out of sample forecasts for the desired variable and the desired periods.}
\usage{
loop_over_forecastDFM(X, variable, r, p, fcst_for_periods)
}
\arguments{
\item{X}{A matrix or data frame}

\item{variable}{The variable to be forecast}

\item{r}{The number of factors to be estimated}

\item{p}{The lag length of the ADL model}

\item{fcst_for_periods}{The periods to be forecast}
}
\value{
The one step ahead forecasts and the real values of the supplied variable
}
\description{
A function to estimate a Dynamic Factor Model and use it to perform pseudo out of sample forecasts for the desired variable and the desired periods.
}
