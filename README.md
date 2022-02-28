# RootId
R code, scripts, and data for unit root identification methodology needed to replicate the empirical results of the paper "Identification of the Differencing Operator of a Non-stationary Time Series via Testing for Zeroes in the Spectral Density".
 
 ## Contents
 1. README.txt: contains overview of how to run the analyses.
 2. ARMAauto.r: R function to compute autocovariances.
 3. flipIt.r: R function to flip explosive autoregressive roots.
 4. Functions.r: other auxiliary R functions.
 5. mvar.loadcon.r: R function to do data cleaning for construction series.
 6. mvar.loadret.r: R function to do data cleaning for retail series. 
 7. root.sigtest.ci.r: R function to apply Root test of "A diagnostic for seasonality based upon polynomial  roots of ARMA models", Journal of Official Statistics 2021.
 8. seas.proxim.r: R function needed to apply Root test.
 9. DataAnalysisTwoSeries.r: R script for the empirical analysis.
 10. Series1Series2.txt: data file in 2-column format.
