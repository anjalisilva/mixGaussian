#' Clustering Using Mixtures of Multivariate Gaussian Via EM
#'
#' Performs clustering using mixtures of multivariate Gaussian
#' distribution with expectation-maximization (EM) for parameter
#' estimation. Model selection is performed using AIC, AIC3,
#' BIC and ICL.
#'
#' @param dataset A dataset of class matrix such that rows correspond to
#'    observations and columns correspond to variables. The dataset have
#'    dimensions n x d, where n is the total number of observations and d
#'    is the dimensionality. If rowSums are zero, these rows will be removed
#'    prior to cluster analysis.
#' @param membership A numeric vector of length nrow(dataset) containing the
#'    cluster membership of each observation. If not available, leave as "none".
#' @param gmin A positive integer specifying the minimum number of components
#'    to be considered in the clustering run.
#' @param gmax A positive integer, >= gmin, specifying the maximum number of
#'    components to be considered in the clustering run.
#' @param initMethod An algorithm for initialization. Current options are
#'    "kmeans", "random", "medoids", "clara", or "fanny". Default is "kmeans".
#' @param nInitIterations A positive integer or zero, specifying the number
#'    of initialization runs to be performed. This many runs, each with 10
#'    iterations, will be performed via mixGaussianClust and values from
#'    the run with highest log-likelihood will be used as initialization
#'    values. Default is 2.
#'
#' @return Returns an S3 object of class mplnVariational with results.
#' \itemize{
#'   \item dataset - The input dataset on which clustering is performed.
#'   \item dimensionality - Dimensionality of the input dataset.
#'   \item gmin - Minimum number of components/clusters considered in the clustering
#'      run.
#'   \item gmax - Maximum number of components/clusters considered in the clustering
#'      run.
#'   \item initalizationMethod - Method used for initialization.
#'   \item allResults - A list with all results.
#'   \item logLikelihood - A vector with value of final log-likelihoods for
#'      each component/cluster size.
#'   \item numbParameters - A vector with number of parameters for each
#'      component/cluster size.
#'   \item trueLabels - The vector of true labels, if provided by user.
#'   \item ICLresults - A list with all ICL model selection results.
#'   \item BICresults - A list with all BIC model selection results.
#'   \item AICresults - A list with all AIC model selection results.
#'   \item AIC3results - A list with all AIC3 model selection results.
#'   \item totalTime - Total time used for clustering and model selection.
#' }
#'
#' @examples
#' # Generating simulated data
#' G <- 2 # number of true clusters/components
#' dimension <- 6
#' nObservations <- 100
#' piGTrue <- c(0.8, 0.2)
#'
#' set.seed(1234)
#' mean1 <- rep(1, dimension)
#' mean2 <- rep(4, dimension)
#' sigma1 <- diag(dimension) * 2
#' sigma2 <- diag(dimension) * 2
#'
#' # library(mvtnorm)
#' component1 <- mvtnorm::rmvnorm(nObservations * piGTrue[1], mean = mean1, sigma = sigma1)
#' dim(component1)
#' component2 <- mvtnorm::rmvnorm(nObservations * piGTrue[2], mean = mean2, sigma = sigma2)
#' dim(component2)
#' dataset <- rbind(component1, component2)
#' dim(dataset) # 100   6
#'
#' # Visualize data
#' # pairs(dataset, col = c(rep(2, nObservations * piGTrue[1]), rep(3, nObservations * piGTrue[2])))
#'
#' # Cluster data
#' clustOutput <- mixGaussian::mixGaussianEM(dataset = dataset,
#'                                           membership = "none",
#'                                           gmin = 1,
#'                                           gmax = 5,
#'                                           initMethod = "kmeans",
#'                                           nInitIterations = 1)
#'
#' @author {Anjali Silva, \email{anjali.silva@uhnresearch.ca}}
#'
#' @references
#' Dempster, A. P., N. M. Laird, and D. B. Rubin (1977). Maximum likelihood from incomplete
#' data via the EM algorithm. \emph{Journal of the Royal Statistical Society: Series B} 39, 1–38.
#'
#' @export
#' @import stats
#' @importFrom mclust map
#' @importFrom mclust unmap
#' @import mvtnorm
#' @import cluster
#'
mixGaussianEM <- function(dataset,
                          membership = "none",
                          gmin,
                          gmax,
                          initMethod = "random",
                          nInitIterations = 2) {

  initialTime <- proc.time()

  # Performing checks
  if (typeof(dataset) != "double") {
    stop("Dataset type needs to be double")
  }

  if (is.matrix(dataset) != TRUE) {
    stop("Dataset needs to be a matrix.")
  }

  if (any(colSums(dataset) <= 0)) {
    stop("Column sums cannot be less than or equal to 0. Double check dataset.")
  }

  dimensionality <- ncol(dataset)
  nObservations <- nrow(dataset)

  if(is.numeric(gmin) != TRUE || is.numeric(gmax) != TRUE) {
    stop("Class of gmin and gmin should be numeric.")
  }

  if (gmax < gmin) {
    stop("gmax cannot be less than gmin.")
  }

  if (gmax > nObservations) {
    stop("gmax cannot be larger than nrow(dataset).")
  }

  if(all(membership != "none") && is.numeric(membership) != TRUE) {
    stop("membership should be a numeric vector containing the
      cluster membership. Otherwise, leave as 'none'.")
  }

  if(all(membership != "none") &&
     all((diff(sort(unique(membership))) == 1) != TRUE) ) {
    stop("Cluster memberships in the membership vector
      are missing a cluster, e.g. 1, 3, 4, 5, 6 is missing cluster 2.")
  }

  if(all(membership != "none") && length(membership) != nObservations) {
    stop("membership should be a numeric vector, where length(membership)
      should equal the number of observations. Otherwise, leave as 'none'.")
  }

  # Remove rows with only zeros, if present
  removezeros <- removeZeroCounts(dataset = dataset, membership = membership)
  dataset <- removezeros$dataset
  membership <- removezeros$membership
  dimensionality <- ncol(dataset)
  nObservations <- nrow(dataset)

  if (is.character(initMethod) == TRUE) {
    initMethodsUsed <- c("kmeans", "random", "medoids", "clara", "fanny")
    if(all((initMethod == initMethodsUsed) == FALSE)) {
      stop("initMethod should of class character, specifying
        either: kmeans, random, medoids, clara, or fanny.")
    }
  } else if (is.character(initMethod) != TRUE) {
    stop("initMethod should of class character, specifying
        either: kmeans, random, medoids, clara, or fanny.")
  }

  if (is.numeric(nInitIterations) != TRUE) {
    stop("nInitIterations should be positive integer or zero, specifying
      the number of initialization runs to be considered.")
  }




  clusterResults <- list() # to save cluster output
  for (gmodel in seq_along(1:(gmax - gmin + 1))) {

    if(length(1:(gmax - gmin + 1)) == gmax) {
      clustersize <- gmodel
    } else if(length(1:(gmax - gmin + 1)) < gmax) {
      clustersize <- seq(gmin, gmax, 1)[gmodel]
    }
    cat("\n Running for g =", clustersize)
    clusterResults[[gmodel]] <- mixGaussianClust(dataset = dataset,
                                                 initMethod = initMethod,
                                                 nInitIterations = nInitIterations,
                                                 G = clustersize,
                                                 maxIterations = 1000)
  }

  names(clusterResults) <- paste0(rep("G=", length(seq(gmin, gmax, 1))), seq(gmin, gmax, 1))

  BIC <- ICL <- AIC <- AIC3 <- Djump <- DDSE <- nParameters <- logLikelihood <- vector()

  for(g in seq_along(1:(gmax - gmin + 1))) {
    # save the final log-likelihood

    if(length(1:(gmax - gmin + 1)) == gmax) {
      clustersize <- g
    } else if(length(1:(gmax - gmin + 1)) < gmax) {
      clustersize <- seq(gmin, gmax, 1)[g]
    }

    logLikelihood[g] <- unlist(utils::tail(clusterResults[[g]]$logLikelihood, n = 1))

    nParameters[g] <- calcParameters(numberG = clustersize, dimensionality = dimensionality)

    if (g == max(1:(gmax - gmin + 1))) { # starting model selection
      bic <- BICFunction(logLikelihood = logLikelihood,
                         nParameters = nParameters,
                         nObservations = nObservations,
                         clusterRunOutput = clusterResults,
                         gmin = gmin,
                         gmax = gmax)

      icl <- ICLFunction(logLikelihood = logLikelihood,
                         nParameters = nParameters,
                         nObservations = nObservations,
                         gmin = gmin,
                         gmax = gmax,
                         clusterRunOutput = clusterResults)

      aic <- AICFunction(logLikelihood = logLikelihood,
                         nParameters = nParameters,
                         clusterRunOutput = clusterResults,
                         gmin = gmin,
                         gmax = gmax)

      aic3 <- AIC3Function(logLikelihood = logLikelihood,
                           nParameters = nParameters,
                           clusterRunOutput = clusterResults,
                           gmin = gmin,
                           gmax = gmax)
    }

  }

  finalTime <- proc.time() - initialTime

  RESULTS <- list(dataset = dataset,
                  dimensionality = dimensionality,
                  gmin = gmin,
                  gmax = gmax,
                  initalizationMethod = initMethod,
                  allResults = clusterResults,
                  logLikelihood = logLikelihood,
                  numbParameters = nParameters,
                  trueLabels = membership,
                  ICLresults = icl,
                  BICresults = bic,
                  AICresults = aic,
                  AIC3results = aic3,
                  totalTime = finalTime)

  class(RESULTS) <- "mixGaussianEM"
  return(RESULTS)

}



mixGaussianClust <- function(dataset,
                             G,
                             initMethod,
                             nInitIterations,
                             maxIterations) {

  dimensionality <- ncol(dataset)
  nObservations <- nrow(dataset)

  # If no intialization is requested by user
  if (nInitIterations == 0) {

    # basic initialization performed for parameters
    mu <- sigma <- num <- list()

    # kmeans initialization
    zValue <- mclust::unmap(kmeans(x = dataset, centers = G, nstart = 100)$cluster)

    # if z generated doesn't add up to nObservations, then use random initialization
    if (sum(colSums(zValue)) != nObservations) {
      zValue <- t(stats::rmultinom(nObservations, size = 1,
                                   prob = rep(1 / G, G)))
    }
    # if z generated has less columns than numbG, then use random initialization
    if(ncol(zValue) < G) {
      zValue <- t(stats::rmultinom(nObservations, size = 1,
                                   prob = rep(1 / G, G)))
    }

    piG <- colSums(zValue) / nObservations


    for (g in 1:G) {
      obs <- which(zValue[ , g] == 1)
      mu[[g]] <- colMeans(dataset[obs, ]) # starting value for mu
      sigma[[g]] <- var(dataset[obs, ]) # starting value for sample covariance matrix
    }

  } else if (nInitIterations != 0) {
    # if initialization is requested by user

    initializationResults <- mixGaussianInit(dataset = dataset,
                                             numbG = G,
                                             initMethod = initMethod,
                                             nInitIterations = nInitIterations)

    mu <- initializationResults$mu
    sigma <- initializationResults$sigma
    zValue <- initializationResults$zValue
    piG <- colSums(zValue) / nObservations

    # other variables
    num <- list()
  }



  # Start clustering
  itOuter <- 1
  aloglik <- logLikelihood <- NULL
  conv <- aloglik[c(1, 2, 3)] <- 0

  while(! conv) {

    # Updating mu
    for(g in 1:G){
      mu[[g]] <- colSums(zValue[ , g] * dataset) / sum(zValue[ , g])
    }

    # Updating covariance matrix
    for(g in 1:G){
      for(i in 1:nObservations){
        num[[i]] <- zValue[i, g]*((dataset[i, ] - mu[[g]])%*%(t(dataset[i, ] - mu[[g]])))
      }
      sigma[[g]] <- Reduce('+', num)/(colSums(zValue)[g])
    }

    # Update pig
    piG <- colSums(zValue) / nObservations


    # Update zvalue
    tempZ <- matrix(0, nObservations, G) # Temporary group membership (z) values
    for(g in 1:G) {
      tempZ[, g] <- piG[g] * mvtnorm::dmvnorm(x = dataset,
                                             mean = mu[[g]],
                                             sigma = sigma[[g]])
    }

    # Calculate zValue value
    # check which tempZ == 0 and rowSums(tempZ)==0 and which of these
    # have both equalling to 0 (because 0/0 = NaN)
    if (G == 1) {
      errorpossible <- Reduce(intersect,
                              list(which(tempZ == 0),
                                   which(rowSums(tempZ) == 0)))
      tempZ[errorpossible] <- 1e-100
      zValue <- tempZ / rowSums(tempZ)
    } else {

      # check for error, if rowsums are zero
      rowSumsZero <- which(rowSums(tempZ) == 0)
      if(length(rowSumsZero) > 1) {
        tempZ[rowSumsZero, ] <- mclust::unmap(stats::kmeans(log(dataset + 1 / 6),
                                                           centers = G,
                                                           nstart = 100)$cluster)[rowSumsZero, ]
        zValue <- tempZ / rowSums(tempZ)
      } else {
        zValue <- tempZ / rowSums(tempZ)
      }
    }



    # Calculate log-likelihood
    logLikelihood[itOuter] <- sum(log(rowSums(tempZ)))


    # Stopping criterion
    if (itOuter > 2) {
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

  # Naming parameters
  names(mu) <- names(sigma) <- paste0(rep("G=", G), 1:G)

  # Saving results for output
  Results <- list(piG = piG,
                  mu = mu,
                  sigma = sigma,
                  probaPost = zValue,
                  clusterlabels = mclust::map(zValue),
                  logLikelihood = logLikelihood)

  class(Results) <- "mixGaussianClustering"
  return(Results)
}


mixGaussianInit <- function(dataset,
                            numbG,
                            initMethod,
                            nInitIterations) {

  dimensionality <- ncol(dataset)
  nObservations <- nrow(dataset)

  zValue <- initRuns <- list()
  logLinit <- vector()

  # Internal function for random initialization
  randomInitfunction <- function(numbG, nObservations) {
    if(numbG == 1) { # generating zValue if g = 1
      zValue <- as.matrix(rep.int(1, times = nObservations),
                          ncol = numbG,
                          nrow = nObservations)
    } else { # generating zValue if g > 1
      zValueConv <- 0
      while(! zValueConv) {
        # ensure that dimension of zValue is same as G (i.e.,
        # if one column contains all 0s, then generate zValue again)
        zValue <- t(stats::rmultinom(nObservations, size = 1,
                                     prob = rep(1 / numbG, numbG)))
        if(length(which(colSums(zValue) > 0)) == numbG) {
          zValueConv <- 1
        }
      }
    }
    return(zValue)
  }

  for(iterations in seq_along(1:nInitIterations)) {
    # setting seed, to ensure if multiple iterations are selected by
    # user, then each run will give a different result.
    set.seed(iterations)
    if (initMethod == "kmeans" | is.na(initMethod)) {
      zValue[[iterations]] <- mclust::unmap(stats::kmeans(dataset,
                                                          centers = numbG, nstart = 100)$cluster)
      # if z generated has less columns than numbG, then use random initialization
      if(ncol(zValue[[iterations]]) < numbG) {
        zValue[[iterations]] <- randomInitfunction(numbG = numbG,
                                                   nObservations = nObservations)
      }

    } else if (initMethod == "random") {
      zValue[[iterations]] <- randomInitfunction(numbG = numbG,
                                                 nObservations = nObservations)

    } else if (initMethod == "medoids") {
      zValue[[iterations]] <- mclust::unmap(cluster::pam(dataset,
                                                         k = numbG,  cluster.only = TRUE))
      # if z generated has less columns than numbG, then use random initialization
      if(ncol(zValue[[iterations]]) < numbG) {
        zValue[[iterations]] <- randomInitfunction(numbG = numbG,
                                                   nObservations = nObservations)
      }

    } else if (initMethod == "clara") {
      zValue[[iterations]] <- mclust::unmap(cluster::clara(dataset,
                                                           k = numbG)$cluster)
      # if z generated has less columns than numbG, then use random initialization
      if(ncol(zValue[[iterations]]) < numbG) {
        zValue[[iterations]] <- randomInitfunction(numbG = numbG, nObservations = nObservations)
      }

    } else if (initMethod == "fanny") {
      zValue[[iterations]] <- mclust::unmap(cluster::fanny(dataset,
                                                           k = numbG,
                                                           memb.exp = numbG,
                                                           cluster.only = TRUE)$clustering)
      # if z generated has less columns than numbG, then use random initialization
      if(ncol(zValue[[iterations]]) < numbG) {
        zValue[[iterations]] <- randomInitfunction(numbG = numbG,
                                                   nObservations = nObservations)
      }
    }

    # maxIterations set to 10 for initialization
    initRuns[[iterations]] <- varMPLNInitClust(dataset = dataset,
                                               G = numbG,
                                               zValue = zValue[[iterations]],
                                               maxIterations = 10)

    logLinit[iterations] <-
      unlist(utils::tail((initRuns[[iterations]]$logLikelihood), n = 1))
  }

  # select the initialization run with highest loglikelihood
  initializationResults <- initRuns[[which(logLinit == max(logLinit, na.rm = TRUE))[1]]]

  class(initializationResults) <- "mixGaussianInit"
  return(initializationResults)
}


varMPLNInitClust <- function(dataset,
                             G,
                             zValue,
                             normFactors,
                             maxIterations = 10) {

  dimensionality <- ncol(dataset)
  nObservations <- nrow(dataset)

  # basic initialization performed for parameters
  mu <- sigma <- num <- list()

  # kmeans initialization
  zValue <- mclust::unmap(kmeans(x = dataset, centers = G, nstart = 100)$cluster)

  # if z generated doesn't add up to nObservations, then use random initialization
  if (sum(colSums(zValue)) != nObservations) {
    zValue <- t(stats::rmultinom(nObservations, size = 1,
                                 prob = rep(1 / G, G)))
  }
  # if z generated has less columns than numbG, then use random initialization
  if(ncol(zValue) < G) {
    zValue <- t(stats::rmultinom(nObservations, size = 1,
                                 prob = rep(1 / G, G)))
  }

  piG <- colSums(zValue) / nObservations


  for (g in 1:G) {
    obs <- which(zValue[ , g] == 1)
    mu[[g]] <- colMeans(dataset[obs, ]) # starting value for mu
    sigma[[g]] <- var(dataset[obs, ]) # starting value for sample covariance matrix
  }


  # Start clustering
  itOuter <- 1
  aloglik <- logLikelihood <- NULL
  conv <- aloglik[c(1, 2, 3)] <- 0

  while(! conv) {

    # Updating mu
    for(g in 1:G){
      mu[[g]] <- colSums(zValue[ , g] * dataset) / sum(zValue[ , g])
    }

    # Updating covariance matrix
    for(g in 1:G){
      for(i in 1:nObservations){
        num[[i]] <- zValue[i, g]*((dataset[i, ] - mu[[g]])%*%(t(dataset[i, ] - mu[[g]])))
      }
      sigma[[g]] <- Reduce('+', num)/(colSums(zValue)[g])
    }

    # Update pig
    piG <- colSums(zValue) / nObservations


    # Update zvalue
    tempZ <- matrix(0, nObservations, G) # Temporary group membership (z) values
    for(g in 1:G) {
      tempZ[, g] <- piG[g] * mvtnorm::dmvnorm(x = dataset,
                                              mean = mu[[g]],
                                              sigma = sigma[[g]])
    }

    # Calculate zValue value
    # check which tempZ == 0 and rowSums(tempZ)==0 and which of these
    # have both equalling to 0 (because 0/0 = NaN)
    if (G == 1) {
      errorpossible <- Reduce(intersect,
                              list(which(tempZ == 0),
                                   which(rowSums(tempZ) == 0)))
      tempZ[errorpossible] <- 1e-100
      zValue <- tempZ / rowSums(tempZ)
    } else {

      # check for error, if rowsums are zero
      rowSumsZero <- which(rowSums(tempZ) == 0)
      if(length(rowSumsZero) > 1) {
        tempZ[rowSumsZero, ] <- mclust::unmap(stats::kmeans(log(dataset + 1 / 6),
                                                            centers = G,
                                                            nstart = 100)$cluster)[rowSumsZero, ]
        zValue <- tempZ / rowSums(tempZ)
      } else {
        zValue <- tempZ / rowSums(tempZ)
      }
    }



    # Calculate log-likelihood
    logLikelihood[itOuter] <- sum(log(rowSums(tempZ)))


    # Stopping criterion
    if (itOuter > 2) {
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

  # Naming parameters
  names(mu) <- names(sigma) <- paste0(rep("G=", G), 1:G)

  # Saving results for output
  Results <- list(piG = piG,
                  mu = mu,
                  sigma = sigma,
                  zValue = zValue,
                  clusterlabels = mclust::map(zValue),
                  logLikelihood = logLikelihood)

  class(Results) <- "varMPLNInitClust"
  return(Results)
}


# Calculate and remove rows with zeros
#' @author {Anjali Silva, \email{anjali.silva@uhnresearch.ca}}
removeZeroCounts <- function(dataset,
                             membership = "none") {

  zeroSUMrows <- which(rowSums(dataset) == 0)

  if (length(zeroSUMrows) > 0 && is.numeric(membership) == TRUE) {
    dataset <- dataset[- zeroSUMrows, ]
    membership <- membership[- zeroSUMrows]
  } else if(length(zeroSUMrows) > 0 && all(membership == "none")) {
    dataset <- dataset[- zeroSUMrows, ]
    membership <- "none"
  }

  RESULTS <- list(dataset = dataset,
                  membership = membership)
  class(RESULTS) <- "ZerosRemoved"
  return(RESULTS)
}


# Parameter calculation
calcParameters <- function(numberG,
                           dimensionality) {

  muPara <- dimensionality * numberG

  sigmaPara <- (dimensionality * ((dimensionality + 1) / 2)) * numberG
  # because if you have numberG-1 parameters,
  # you can do 1-these to get the last one

  piPara <- numberG - 1

  # total parameters
  paraTotal <- muPara + sigmaPara + piPara

  return(paraTotal)
  # Developed by Anjali Silva
}


#' Model Selection Via Akaike Information Criterion
#'
#' Performs model selection using Akaike Information Criterion (AIC).
#' Formula: - 2 * logLikelihood + 2 * nParameters.
#'
#' @param logLikelihood A vector with value of final log-likelihoods for
#'      each cluster size.
#' @param nParameters A vector with number of parameters for each
#'      cluster size.
#' @param clusterRunOutput Output from mixGaussianEM, if available. Default
#'    value is NA. If provided, the vector of cluster labels obtained by
#'    mclust::map() for best model will be provided in the output.
#' @param gmin A positive integer specifying the minimum number of components
#'    to be considered in the clustering run.
#' @param gmax A positive integer, >gmin, specifying the maximum number of
#'    components to be considered in the clustering run.
#'
#' @return Returns an S3 object of class MPLN with results.
#' \itemize{
#'   \item allAICvalues - A vector of AIC values for each cluster size.
#'   \item AICmodelselected - An integer specifying model selected by AIC.
#'   \item AICmodelSelectedLabels - A vector of integers specifying cluster labels
#'     for the model selected. Only provided if user input clusterRunOutput.
#'   \item AICMessage - A character vector indicating if spurious clusters are
#'     detected. Otherwise, NA.
#' }
#'
#' @examples
#'
#' # Generating simulated data
#' G <- 2 # number of true clusters/components
#' dimension <- 6
#' nObservations <- 100
#' piGTrue <- c(0.8, 0.2)
#'
#' set.seed(1234)
#' mean1 <- rep(1, dimension)
#' mean2 <- rep(4, dimension)
#' sigma1 <- diag(dimension) * 0.2
#' sigma2 <- diag(dimension) * 0.5
#'
#' # library(mvtnorm)
#' component1 <- mvtnorm::rmvnorm(nObservations * piGTrue[1], mean = mean1, sigma = sigma1)
#' component2 <- mvtnorm::rmvnorm(nObservations * piGTrue[2], mean = mean2, sigma = sigma2)
#' dataset <- rbind(component1, component2)
#'
#' # Cluster data
#' clustOutput <- mixGaussian::mixGaussianEM(dataset = dataset,
#'                                           membership = "none",
#'                                           gmin = 1,
#'                                           gmax = 2,
#'                                           initMethod = "kmeans",
#'                                           nInitIterations = 1)
#'
#' # Model selection
#' AICmodel <- mixGaussian::AICFunction(logLikelihood = clustOutput$logLikelihood,
#'                                      nParameters = clustOutput$numbParameters,
#'                                      clusterRunOutput = clustOutput$allResults,
#'                                      gmin = clustOutput$gmin,
#'                                      gmax = clustOutput$gmax)
#'
#' @author {Anjali Silva, \email{anjali.silva@uhnresearch.ca}}
#'
#' @references
#'
#' Akaike, H. (1973). Information theory and an extension of the maximum likelihood
#' principle. In \emph{Second International Symposium on Information Theory}, New York, NY,
#' USA, pp. 267–281. Springer Verlag.
#'
#' @export
#'
AICFunction <- function(logLikelihood,
                        nParameters,
                        clusterRunOutput = NA,
                        gmin,
                        gmax) {

  # Performing checks
  if(is.numeric(gmin) != TRUE || is.numeric(gmax) != TRUE) {
    stop("Class of gmin and gmin should be numeric.")
  }

  if (gmax < gmin) {
    stop("gmax cannot be less than gmin.")
  }

  if(is.numeric(nParameters) != TRUE) {
    stop("nParameters should be a vector of integers indicating
      number of parameters for each cluster size.")
  }

  if(is.numeric(logLikelihood) != TRUE) {
    stop("logLikelihood should be a vector of numeric values.")
  }

  if(length(logLikelihood) != (gmax - gmin + 1)) {
    stop("logLikelihood should be a vector of length (gmax - gmin + 1).")
  }

  AIC <- - 2 * logLikelihood + 2 * nParameters
  AICmodel <- seq(gmin, gmax, 1)[grep(min(AIC, na.rm = TRUE), AIC)]
  AICMessage <- NA # For spurious clusters


  # Obtain cluster labels for best model if clusterRunOutput is provided
  if(all(is.na(clusterRunOutput)) == FALSE) {
      AICmodelLabels <- clusterRunOutput[[grep(min(AIC, na.rm = TRUE), AIC)]]$clusterlabels

     # Check for spurious clusters, only possible if cluster labels provided
     if (max(AICmodelLabels) != AICmodel) {
      AICmodel <- max(AICmodelLabels)
      AICMessage <- "Spurious or empty cluster resulted."
    }
  } else {
    AICmodelLabels <- "clusterRunOutput not provided"
  }

  AICresults<-list(allAICvalues = AIC,
                   AICmodelselected = AICmodel,
                   AICmodelSelectedLabels = AICmodelLabels,
                   AICMessage = AICMessage)
  class(AICresults) <- "AIC"
  return(AICresults)
}


#' Model Selection Via Akaike Information Criterion by Bozdogan (1994)
#'
#' Performs model selection using Akaike Information Criterion by
#' Bozdogan (1994), called AIC3. Formula: - 2 * logLikelihood + 3 * nParameters.
#'
#' @param logLikelihood A vector with value of final log-likelihoods for
#'      each cluster size.
#' @param nParameters A vector with number of parameters for each
#'      cluster size.
#' @param clusterRunOutput Output from mixGaussianEM, if available. Default
#'    value is NA. If provided, the vector of cluster labels obtained by
#'    mclust::map() for best model will be provided in the output.
#' @param gmin A positive integer specifying the minimum number of components
#'    to be considered in the clustering run.
#' @param gmax A positive integer, >gmin, specifying the maximum number of
#'    components to be considered in the clustering run.
#'
#' @return Returns an S3 object of class MPLN with results.
#' \itemize{
#'   \item allAIC3values - A vector of AIC3 values for each cluster size.
#'   \item AIC3modelselected - An integer specifying model selected by AIC3.
#'   \item AIC3modelSelectedLabels - A vector of integers specifying cluster labels
#'     for the model selected. Only provided if user input clusterRunOutput.
#'   \item AIC3Message - A character vector indicating if spurious clusters are
#'     detected. Otherwise, NA.
#' }
#'
#' @examples
#'
#' # Generating simulated data
#' G <- 2 # number of true clusters/components
#' dimension <- 6
#' nObservations <- 100
#' piGTrue <- c(0.8, 0.2)
#'
#' set.seed(1234)
#' mean1 <- rep(1, dimension)
#' mean2 <- rep(4, dimension)
#' sigma1 <- diag(dimension) * 0.2
#' sigma2 <- diag(dimension) * 0.5
#'
#' # library(mvtnorm)
#' component1 <- mvtnorm::rmvnorm(nObservations * piGTrue[1], mean = mean1, sigma = sigma1)
#' component2 <- mvtnorm::rmvnorm(nObservations * piGTrue[2], mean = mean2, sigma = sigma2)
#' dataset <- rbind(component1, component2)
#'
#' # Cluster data
#' clustOutput <- mixGaussian::mixGaussianEM(dataset = dataset,
#'                                           membership = "none",
#'                                           gmin = 1,
#'                                           gmax = 2,
#'                                           initMethod = "kmeans",
#'                                           nInitIterations = 1)
#'
#' # Model selection
#' AIC3model <- mixGaussian::AIC3Function(logLikelihood = clustOutput$logLikelihood,
#'                                        nParameters = clustOutput$numbParameters,
#'                                        clusterRunOutput = clustOutput$allResults,
#'                                        gmin = clustOutput$gmin,
#'                                        gmax = clustOutput$gmax)
#'
#' @author {Anjali Silva, \email{anjali.silva@uhnresearch.ca}}
#'
#' @references
#'
#' Akaike, H. (1973). Information theory and an extension of the maximum likelihood
#' principle. In \emph{Second International Symposium on Information Theory}, New York, NY,
#' USA, pp. 267–281. Springer Verlag.
#'
#' #' Bozdogan, H. (1994). Mixture-model cluster analysis using model selection criteria
#' and a new informational measure of complexity. In \emph{Proceedings of the First US/Japan
#' Conference on the Frontiers of Statistical Modeling: An Informational Approach:
#' Volume 2 Multivariate Statistical Modeling}, pp. 69–113. Dordrecht: Springer Netherlands.
#'
#' @export
#'
AIC3Function <- function(logLikelihood,
                         nParameters,
                         clusterRunOutput = NA,
                         gmin,
                         gmax) {

  # Performing checks
  if(is.numeric(gmin) != TRUE || is.numeric(gmax) != TRUE) {
    stop("Class of gmin and gmin should be numeric.")
  }

  if (gmax < gmin) {
    stop("gmax cannot be less than gmin.")
  }

  if(is.numeric(nParameters) != TRUE) {
    stop("nParameters should be a vector of integers indicating
      number of parameters for each cluster size.")
  }

  if(is.numeric(logLikelihood) != TRUE) {
    stop("logLikelihood should be a vector of numeric values.")
  }

  if(length(logLikelihood) != (gmax - gmin + 1)) {
    stop("logLikelihood should be a vector of length (gmax - gmin + 1).")
  }

  AIC3 <- - 2 * logLikelihood + 3 * nParameters
  AIC3model <- seq(gmin, gmax, 1)[grep(min(AIC3,na.rm = TRUE), AIC3)]
  AIC3Message <- NA # For spurious clusters

  # Obtain cluster labels for best model if clusterRunOutput is provided
  if(all(is.na(clusterRunOutput)) == FALSE) {
    AIC3modelLabels <- clusterRunOutput[[grep(min(AIC3,na.rm = TRUE), AIC3)]]$clusterlabels

    # Check for spurious clusters, only possible if cluster labels provided
    if (max(AIC3modelLabels) != AIC3model) {
      AIC3model <- max(AIC3modelLabels)
      AIC3Message <- "Spurious or empty cluster resulted."
    }
  } else {
    AIC3modelLabels <- "clusterRunOutput not provided"
  }

  AIC3results <- list(allAIC3values = AIC3,
                      AIC3modelselected = AIC3model,
                      AIC3modelSelectedLabels = AIC3modelLabels,
                      AIC3Message = AIC3Message)
  class(AIC3results) <- "AIC3"
  return(AIC3results)
}



#' Model Selection Via Bayesian Information Criterion
#'
#' Performs model selection using Bayesian Information Criterion (BIC) by
#' Schwarz (1978). Formula: - 2 * logLikelihood + (nParameters * log(nObservations)).
#'
#' @param logLikelihood A vector with value of final log-likelihoods for
#'      each cluster size.
#' @param nParameters A vector with number of parameters for each
#'      cluster size.
#' @param nObservations A positive integer specifying the number of observations
#'      in the dataset analyzed.
#' @param clusterRunOutput Output from mixGaussianEM, if available. Default value
#'    is NA. If provided, the vector of cluster labels obtained by mclust::map()
#'    for best model will be provided in the output.
#' @param gmin A positive integer specifying the minimum number of components
#'    to be considered in the clustering run.
#' @param gmax A positive integer, >gmin, specifying the maximum number of
#'    components to be considered in the clustering run.
#'
#' @return Returns an S3 object of class MPLN with results.
#' \itemize{
#'   \item allBICvalues - A vector of BIC values for each cluster size.
#'   \item BICmodelselected - An integer specifying model selected by BIC
#'   \item BICmodelSelectedLabels - A vector of integers specifying cluster labels
#'     for the model selected. Only provided if user input clusterRunOutput.
#'   \item BICMessage - A character vector indicating if spurious clusters are
#'     detected. Otherwise, NA.
#' }
#'
#' @examples
#'
#' # Generating simulated data
#' G <- 2 # number of true clusters/components
#' dimension <- 6
#' nObservations <- 100
#' piGTrue <- c(0.8, 0.2)
#'
#' set.seed(1234)
#' mean1 <- rep(1, dimension)
#' mean2 <- rep(4, dimension)
#' sigma1 <- diag(dimension) * 0.2
#' sigma2 <- diag(dimension) * 0.5
#'
#' # library(mvtnorm)
#' component1 <- mvtnorm::rmvnorm(nObservations * piGTrue[1], mean = mean1, sigma = sigma1)
#' component2 <- mvtnorm::rmvnorm(nObservations * piGTrue[2], mean = mean2, sigma = sigma2)
#' dataset <- rbind(component1, component2)
#'
#' # Cluster data
#' clustOutput <- mixGaussian::mixGaussianEM(dataset = dataset,
#'                                           membership = "none",
#'                                           gmin = 1,
#'                                           gmax = 2,
#'                                           initMethod = "kmeans",
#'                                           nInitIterations = 1)
#'
#' # Model selection
#' BICmodel <- mixGaussian::BICFunction(logLikelihood = clustOutput$logLikelihood,
#'                                      nParameters = clustOutput$numbParameters,
#'                                      nObservations = nrow(clustOutput$dataset),
#'                                      clusterRunOutput = clustOutput$allResults,
#'                                      gmin = clustOutput$gmin,
#'                                      gmax = clustOutput$gmax)
#'
#' @author {Anjali Silva, \email{anjali.silva@uhnresearch.ca}}
#'
#' @references
#'
#' Schwarz, G. (1978). Estimating the dimension of a model. \emph{The Annals of Statistics}
#' 6.
#'
#' @export
#'
BICFunction <- function(logLikelihood,
                        nParameters,
                        nObservations,
                        clusterRunOutput = NA,
                        gmin,
                        gmax) {

  # Performing checks
  if(is.numeric(gmin) != TRUE || is.numeric(gmax) != TRUE) {
    stop("Class of gmin and gmin should be numeric.")
  }

  if (gmax < gmin) {
    stop("gmax cannot be less than gmin.")
  }

  if(is.numeric(nParameters) != TRUE) {
    stop("nParameters should be a vector of integers indicating
      number of parameters for each cluster size.")
  }

  if(is.numeric(logLikelihood) != TRUE) {
    stop("logLikelihood should be a vector of numeric values.")
  }

  if(length(logLikelihood) != (gmax - gmin + 1)) {
    stop("logLikelihood should be a vector of length (gmax - gmin + 1).")
  }

  BIC <- - 2 * logLikelihood + (nParameters * log(nObservations))
  BICmodel <- seq(gmin, gmax, 1)[grep(min(BIC, na.rm = TRUE), BIC)]
  BICMessage <- NA # For spurious clusters

  # Obtain cluster labels for best model if clusterRunOutput is provided
  if(all(is.na(clusterRunOutput)) == FALSE) {
      BICmodelLabels <- clusterRunOutput[[grep(min(BIC, na.rm = TRUE),
                                               BIC)]]$clusterlabels
    # Check for spurious clusters, only possible if cluster labels provided
    if (max(BICmodelLabels) != BICmodel) {
      BICmodel <- max(BICmodelLabels)
      BICMessage <- "Spurious or empty cluster resulted."
    }
  } else {
    BICmodelLabels <- "clusterRunOutput not provided"
  }

  BICresults <- list(allBICvalues = BIC,
                     BICmodelselected = BICmodel,
                     BICmodelSelectedLabels = BICmodelLabels,
                     BICMessage = BICMessage)
  class(BICresults) <- "BIC"
  return(BICresults)
}



#' Model Selection Via Integrated Completed Likelihood
#'
#' Performs model selection using integrated completed likelihood (ICL) by
#' Biernacki et al., (2000).
#'
#' @param logLikelihood A vector with value of final log-likelihoods for
#'    each cluster size.
#' @param nParameters A vector with number of parameters for each
#'    cluster size.
#' @param nObservations A positive integer specifying the number of observations
#'    in the dataset analyzed.
#' @param clusterRunOutput Output from mixGaussianEM function. Either clusterRunOutput
#'    or probaPost must be provided.
#' @param probaPost A list that is length (gmax - gmin + 1) containing posterior
#'    probability at each g, for g = gmin:gmax. This argument is useful if
#'    clustering output have been generated non-serially, e.g., g = 1:5 and
#'    g = 6:10 rather than g = 1:10. Either clusterRunOutput or probaPost
#'    must be provided.
#' @param gmin A positive integer specifying the minimum number of components
#'    to be considered in the clustering run.
#' @param gmax A positive integer, > gmin, specifying the maximum number of
#'    components to be considered in the clustering run.
#'
#' @return Returns an S3 object of class MPLN with results.
#' \itemize{
#'   \item allICLvalues - A vector of ICL values for each cluster size.
#'   \item ICLmodelselected - An integer specifying model selected by ICL.
#'   \item ICLmodelSelectedLabels - A vector of integers specifying cluster labels
#'     for the model selected. Only provided if user input clusterRunOutput.
#'   \item ICLMessage - A character vector indicating if spurious clusters are
#'     detected. Otherwise, NA.
#' }
#'
#' @examples
#'
#' # Generating simulated data
#' G <- 2 # number of true clusters/components
#' dimension <- 6
#' nObservations <- 100
#' piGTrue <- c(0.8, 0.2)
#'
#' set.seed(1234)
#' mean1 <- rep(1, dimension)
#' mean2 <- rep(4, dimension)
#' sigma1 <- diag(dimension) * 0.2
#' sigma2 <- diag(dimension) * 0.5
#'
#' # library(mvtnorm)
#' component1 <- mvtnorm::rmvnorm(nObservations * piGTrue[1], mean = mean1, sigma = sigma1)
#' component2 <- mvtnorm::rmvnorm(nObservations * piGTrue[2], mean = mean2, sigma = sigma2)
#' dataset <- rbind(component1, component2)
#'
#' # Cluster data
#' clustOutput <- mixGaussian::mixGaussianEM(dataset = dataset,
#'                                           membership = "none",
#'                                           gmin = 1,
#'                                           gmax = 2,
#'                                           initMethod = "kmeans",
#'                                           nInitIterations = 1)
#'
#' # Model selection
#' ICLmodel <- mixGaussian::ICLFunction(logLikelihood = clustOutput$logLikelihood,
#'                                      nParameters = clustOutput$numbParameters,
#'                                      nObservations = nrow(clustOutput$dataset),
#'                                      clusterRunOutput = clustOutput$allResults,
#'                                      gmin = clustOutput$gmin,
#'                                      gmax = clustOutput$gmax)
#'
#' @author {Anjali Silva, \email{anjali.silva@uhnresearch.ca}}
#'
#' @references
#'
#' Biernacki, C., G. Celeux, and G. Govaert (2000). Assessing a mixture model for
#' clustering with the integrated classification likelihood. \emph{IEEE Transactions
#' on Pattern Analysis and Machine Intelligence} 22.
#'
#' @export
#' @importFrom mclust unmap
#'
ICLFunction <- function(logLikelihood,
                        nParameters,
                        nObservations,
                        clusterRunOutput = NA,
                        probaPost = NA,
                        gmax,
                        gmin) {

  # Performing checks
  if (gmax < gmin) {
    stop("gmax cannot be less than gmin.")
  }

  if(is.numeric(nParameters) != TRUE) {
    stop("nParameters should be a vector of integers indicating
      number of parameters for each cluster size.")
  }

  if(is.numeric(logLikelihood) != TRUE) {
    stop("logLikelihood should be a vector of numeric values.")
  }

  if(length(logLikelihood) != (gmax - gmin + 1)) {
    stop("logLikelihood should be a vector of length (gmax - gmin + 1).")
  }


  if(all(is.na(probaPost)) != TRUE) {
    if(length(probaPost) != (gmax - gmin + 1)) {
      stop("probaPost must be a list of length (gmax - gmin + 1)
      containing posterior probability at each g.")
    }
  }

  if(all(is.na(clusterRunOutput)) == TRUE && all(is.na(probaPost)) == TRUE) {
    stop("Either clusterRunOutput or probaPost must be provided.")
  }

  BIC <- - 2 * logLikelihood + (nParameters * log(nObservations))

  ICL <- vector()
  forICL <- function(g){sum(log(z[which(mapz[, g] == 1), g]))}

  # if clusterRunOutput is provided by user
  if(all(is.na(clusterRunOutput)) != TRUE) {
    for (g in 1:(gmax - gmin + 1)) {
        z <- clusterRunOutput[[g]]$probaPost
        mapz <- mclust::unmap(clusterRunOutput[[g]]$clusterlabels)
        ICL[g] <- BIC[g] + sum(sapply(1:ncol(mapz), forICL))
    }

    # select best model
    ICLmodel <- seq(gmin, gmax, 1)[grep(min(ICL, na.rm = TRUE), ICL)]

    # select cluster labels for best model
    ICLmodelLabels <- clusterRunOutput[[grep(min(ICL, na.rm = TRUE),
                                               ICL)]]$clusterlabels
  }

  # if probaPost is provided by user
  if(all(is.na(probaPost)) != TRUE) {
    for (g in 1:(gmax - gmin + 1)) {
      z <- probaPost[[g]]
      mapz <- mclust::unmap(mclust::map(probaPost[[g]]))
      forICL <- function(g){sum(log(z[which(mapz[, g] == 1), g]))}
      ICL[g] <- BIC[g] + sum(sapply(1:ncol(mapz), forICL))
    }
    ICLmodel <- seq(gmin, gmax, 1)[grep(min(ICL, na.rm = TRUE), ICL)]
    ICLmodelLabels <- mclust::map(probaPost[[grep(min(ICL, na.rm = TRUE), ICL)]])
  }

  # Check for spurious clusters
  ICLMessage <- NA
  if (max(ICLmodelLabels) != ICLmodel) {
    ICLmodel <- max(ICLmodelLabels)
    ICLMessage <- "Spurious or empty cluster resulted."
  }

  ICLresults <- list(allICLvalues = ICL,
                     ICLmodelselected = ICLmodel,
                     ICLmodelSelectedLabels = ICLmodelLabels,
                     ICLMessage = ICLMessage)
  class(ICLresults) <- "ICL"
  return(ICLresults)
}

# [END]
