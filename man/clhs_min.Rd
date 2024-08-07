% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clhs_min.R
\name{clhs_min}
\alias{clhs_min}
\title{Determine minimum sample size for the clhs algorithm}
\usage{
clhs_min(covs, lq = 0.05, uq = 0.95, vif = FALSE, vif.thresh = 5)
}
\arguments{
\item{covs}{dataframe*, data frame or SpacVector or SpatRaster, environmental covariates to be used for the sample size determination}

\item{lq, }{numeric*, lower quantile of choice for the confidence interval, default value = 0.05}

\item{uq, }{numeric*, upper quantile of choice for the confidence interval, default value = 0.95}

\item{vif, }{logical*, use variance inflation factor to reduce multicolinearity, default = FALSE}

\item{vif.thresh, }{numeric*, threshold to be used for VIF, typically 5 or 10, default = 5}
}
\value{
list, with 7 objects
Object 1. Dataframe 'Sample_Size' returning descriptive statistics for the minimum sample size
Object 2. Dataframe 'Transformations' provides data transformations used to normalize covariate distriburions
Object 3. Dataframe 'Results' provides the minimum sample size for each individual covariate
Object 4. Dataframe 'Bimodality' provides the results of bimodality assessment for covariates deemd to be bimodal
Object 5. List 'Hist_Original' containing histogtram objects for each covariate
Object 6. List 'Hist_Transformed' containing histogram objects for all transformed covariates, including the bimodal covariates whose distributions were split
Object 7. Dataframe with sample plan developed using median of minimum sample size
Object 8. List 'VIF_Results' containing the outputs of the variance inflation factor analysis (optional)
}
\description{
This script is based on the publication:
Determining minimum sample size for the conditioned Latin hypercube sampling algorithm
DOI: https://doi.org/10.1016/j.pedsph.2022.09.001
}
\examples{
data(covs)
out<- clhs_min(covs=covs, lq=0.05,uq=0.95, vif=TRUE, vif.thresh=10)





}
