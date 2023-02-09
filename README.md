# DynamicFactorModel

This is a package that has its origin in my bachelor thesis on "Exploring the Potential of Dynamic Factor Models for High-Dimensional Time Series Forecasting". 

It can be downloaded with

```R
library(devtools)
install_github("linuswolff/DynamicFactorModel")
```

The package has three functions. The first estimates a DFM of **$X$** and forecasts the desired variable, allowing for modification of the number of factors. Criteria for estimating the number of factors as proposed in Bai and Ng (2002) are offered by another. The last performs out of sample forecasts for a desired length of periods. Their documentation can be accessed by using 

```R
?forecastDFM
?ICfactors_mine
?loop_over_forecastDFM
```

The Dynamic Factor Model are estimated by PCA and forecasts are based on the Factor Augmented Auto Regression approach in Stock and Watson (2002).