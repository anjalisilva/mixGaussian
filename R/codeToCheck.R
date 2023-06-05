# Author: Anjali Silva
# Date: June 2023
# Purpose: Develop a discriminant analysis
# program for mixtures of Gaussian. Here, there are
# 100 observations, 90 labelled and 10 unlabeled. The
# 90 labeled observations come from 2 groups.
# The purpose is to assign 10 unlabeled observations to
# 2 existing or more groups.

#### Packages ####
library(mclust)
library(mvtnorm)

#### Example 1: Generating simulated data ####
set.seed(1234)
G <- 2 # number of true clusters/components
dimension <- 6
nObservations <- 100
membershipTrue <- sample(c(1, 2), 100, replace = TRUE)
piGTrue <- table(membershipTrue) / nObservations

# Simulate parameter values
mean1 <- rep(1, dimension)
mean2 <- rep(5, dimension)
sigma1 <- diag(dimension) * 2
sigma2 <- diag(dimension) *1

# library(mvtnorm)
component1 <- mvtnorm::rmvnorm(
                n = as.numeric(nObservations * piGTrue[1]),
                mean = mean1,
                sigma = sigma1)
dim(component1)
component2 <- mvtnorm::rmvnorm(
                n = 100 - as.numeric(nObservations * piGTrue[1]),
                mean = mean2,
                sigma = sigma2)
dim(component2)
dataset <- rbind(component1, component2)
dim(dataset) # 100   6

# Visualize data
pairs(dataset, col = c(rep(2, nObservations * piGTrue[1]),
                       rep(3, nObservations * piGTrue[2])))

# Classification membership vector
membershipClass <- membershipTrue
membershipClass[x = sample(c(1:100),
                 size = 10,
                 replace = FALSE)] <- 0
table(membershipClass)
#  0  1  2
# 10 38 52

#### Begin code ####

dimensionality <- ncol(dataset)
nObservations <- nrow(dataset)
# basic initialization performed for parameters
mu <- sigma <- num <- list()

# initialization
Gtesting <- 2:3 # for uknown memberships, need to classify for 2 group
# or test if there is even a third case
Gknown <- 1:2 # know these exist based on present membership
zValue <- matrix(0, ncol = max(Gtesting), nrow = nObservations)
clsInd <- (membershipClass == 0) # membership unknown index
for (i in 1:nObservations) {
  if(clsInd[i]) { # if membership unknown
    zValue[i, ] <- 1 / max(Gtesting)
  } else {
    zValue[i, membershipClass[i]] <- 1
  }
}
sum(rowSums(zValue)) # check 100 observations
piG <- colSums(zValue) / nObservations
# 0.41333333 0.55333333 0.03333333

for (g in 1:max(Gtesting)) {
  obs <- which(zValue[, g] == 1)
  mu[[g]] <- colMeans(dataset[obs, ]) # starting value for mu
  sigma[[g]] <- var(dataset[obs, ]) # starting value for sample covariance matrix
}


####  Start discriminant analysis ####
itOuter <- 1
aloglik <- logLikelihood <- NULL
conv <- aloglik[c(1, 2, 3)] <- 0
maxIterations <- 1000

while(! conv) {

  # Updating mu using known membership values
  for(g in 1:max(Gtesting)){
    mu[[g]] <- colSums(zValue[ , g] * dataset) / sum(zValue[ , g])
  }

  # Updating covariance matrix
  # using known membership values
  for(g in 1:max(Gtesting)) {
    for(i in 1:nObservations) {
      num[[i]] <- zValue[i, g]*((dataset[i, ] - mu[[g]])%*%(t(dataset[i, ] - mu[[g]])))
    }
    sigma[[g]] <- Reduce('+', num)/(colSums(zValue)[g])
  }



  # Update zvalue for unknown memberships
  # Temporary group membership (z) values

  tempZ <- matrix(0, nObservations, max(Gtesting))
  for(g in 1:max(Gtesting)) {
    tempZ[, g] <- piG[g] * mvtnorm::dmvnorm(x = dataset,
                                            mean = mu[[g]],
                                            sigma = sigma[[g]])
  }

  # Calculate zValue value for unlabeled observations
  zValue <- tempZ / rowSums(tempZ)


  # Calculate log-likelihood
  logLikelihood[itOuter] <- sum(log(rowSums(zValue)))


  # Update pig
  piG <- colSums(zValue) / nObservations

  # Stopping criterion
  if (itOuter > 5) {
    if ((logLikelihood[itOuter - 1] - logLikelihood[itOuter - 2]) == 0) {
      conv <-1
    } else {
      # Aitken stopping criterion
      termAitkens <- (logLikelihood[itOuter]- logLikelihood[itOuter - 1]) /
        (logLikelihood[itOuter - 1] - logLikelihood[itOuter - 2])
      term2Aitkens <- (1 / (1 - termAitkens) * (logLikelihood[itOuter] -
                                                  logLikelihood[itOuter - 1]))
      aloglik[itOuter] <- logLikelihood[itOuter - 1] + term2Aitkens
      if (abs(aloglik[itOuter] - logLikelihood[itOuter - 1]) < 0.001) {
        # If this critera, as per Böhning et al., 1994 is achieved
        # convergence is achieved
        conv <- 1
      } else {
        conv <- conv
      }
    }

  }

  # Update iteration
  itOuter <- itOuter + 1
  if (itOuter == maxIterations) {
    checks <- 1
  }
}


#### Finalize memberships ####
# Assign new memberships to where 0s (no labels) occured
membershipClass[clsInd] <- mclust::map(zValue[clsInd, ])
table(membershipClass)
# 1  2  3
# 39 57  4
Zfinal <- mclust::unmap(membershipClass)
