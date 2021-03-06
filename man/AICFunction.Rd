% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functionsmixGau.R
\name{AICFunction}
\alias{AICFunction}
\title{Model Selection Via Akaike Information Criterion}
\usage{
AICFunction(logLikelihood, nParameters, clusterRunOutput = NA, gmin, gmax)
}
\arguments{
\item{logLikelihood}{A vector with value of final log-likelihoods for
each cluster size.}

\item{nParameters}{A vector with number of parameters for each
cluster size.}

\item{clusterRunOutput}{Output from mixGaussianEM, if available. Default
value is NA. If provided, the vector of cluster labels obtained by
mclust::map() for best model will be provided in the output.}

\item{gmin}{A positive integer specifying the minimum number of components
to be considered in the clustering run.}

\item{gmax}{A positive integer, >gmin, specifying the maximum number of
components to be considered in the clustering run.}
}
\value{
Returns an S3 object of class MPLN with results.
\itemize{
\item allAICvalues - A vector of AIC values for each cluster size.
\item AICmodelselected - An integer specifying model selected by AIC.
\item AICmodelSelectedLabels - A vector of integers specifying cluster labels
for the model selected. Only provided if user input clusterRunOutput.
\item AICMessage - A character vector indicating if spurious clusters are
detected. Otherwise, NA.
}
}
\description{
Performs model selection using Akaike Information Criterion (AIC).
Formula: - 2 * logLikelihood + 2 * nParameters.
}
\examples{

# Generating simulated data
G <- 2 # number of true clusters/components
dimension <- 6
nObservations <- 100
piGTrue <- c(0.8, 0.2)

set.seed(1234)
mean1 <- rep(1, dimension)
mean2 <- rep(4, dimension)
sigma1 <- diag(dimension) * 0.2
sigma2 <- diag(dimension) * 0.5

# library(mvtnorm)
component1 <- mvtnorm::rmvnorm(nObservations * piGTrue[1], mean = mean1, sigma = sigma1)
component2 <- mvtnorm::rmvnorm(nObservations * piGTrue[2], mean = mean2, sigma = sigma2)
dataset <- rbind(component1, component2)

# Cluster data
clustOutput <- mixGaussian::mixGaussianEM(dataset = dataset,
                                          membership = "none",
                                          gmin = 1,
                                          gmax = 2,
                                          initMethod = "kmeans",
                                          nInitIterations = 1)

# Model selection
AICmodel <- mixGaussian::AICFunction(logLikelihood = clustOutput$logLikelihood,
                                     nParameters = clustOutput$numbParameters,
                                     clusterRunOutput = clustOutput$allResults,
                                     gmin = clustOutput$gmin,
                                     gmax = clustOutput$gmax)

}
\references{
Akaike, H. (1973). Information theory and an extension of the maximum likelihood
principle. In \emph{Second International Symposium on Information Theory}, New York, NY,
USA, pp. 267–281. Springer Verlag.
}
\author{
{Anjali Silva, \email{anjali.silva@uhnresearch.ca}}
}
