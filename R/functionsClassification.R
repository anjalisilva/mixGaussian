#' Classification Using Mixtures of Multivariate Gaussian Via EM
#'
#' Performs classification using mixtures of multivariate Gaussian
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
#'    cluster membership of each observation. If not available, leave as "none"
#'    and model-based clustering will be performed. If a vector provided, then
#'    model-based classification is performed. In this latter case, the ith entry
#'    of membership is either zero, indicating that the component membership of
#'    observation i is unknown, or it corresponds to the component membership of
#'    observation i. See Examples below.
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
#' @return Returns an S3 object of class mixGaussianClass with results.
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
#' # Example 1: Generating simulated data
#' set.seed(1234)
#' G <- 2 # number of true clusters/components
#' dimension <- 6
#' nObservations <- 100
#' membershipTrue <- sample(c(1, 2), 100, replace = TRUE)
#' piGTrue <- table(membershipTrue) / nObservations
#'
#' # Simulate parameter values
#' mean1 <- rep(1, dimension)
#' mean2 <- rep(5, dimension)
#' sigma1 <- diag(dimension) * 2
#' sigma2 <- diag(dimension) *1
#'
#' # library(mvtnorm)
#' component1 <- mvtnorm::rmvnorm(
#'                 n = as.numeric(nObservations * piGTrue[1]),
#'                 mean = mean1,
#'                 sigma = sigma1)
#' dim(component1)
#' component2 <- mvtnorm::rmvnorm(
#'                 n = 100 - as.numeric(nObservations * piGTrue[1]),
#'                 mean = mean2,
#'                 sigma = sigma2)
#' dim(component2)
#' dataset <- rbind(component1, component2)
#' dim(dataset) # 100   6
#'
#' # Visualize data
#' pairs(dataset, col = c(rep(2, nObservations * piGTrue[1]),
#'                          rep(3, nObservations * piGTrue[2])))
#'
#' # Classification membership vector
#' membershipClass <- membershipTrue
#' membershipClass[x = sample(c(1:100),
#'                 size = 10,
#'                 replace = FALSE)] <- 0
#' table(membershipClass)
#' #  0  1  2
#' # 10 38 52
#'
#' # Classify data, where cluster membership of observations
#' # with membership indicated by 0 are unknown
#' clustOutput <- mixGaussian::mixGaussianEM(dataset = dataset,
#'                                           membership = membershipClass,
#'                                           gmin = 1,
#'                                           gmax = 5,
#'                                           initMethod = "kmeans",
#'                                           nInitIterations = 1)
#'
#'
#' # Example 2: Access raw data made available with this package
#' # Not run
#' # inputCountsPath <- system.file("extdata", "mixGaussianDataset.csv",
#' #                               package = "mixGaussian")
#' # Read data
#' # exampleData <- read.csv(file = inputCountsPath, header = TRUE)
#' # dim(exampleData) # 1000  6
#' # To see documentation for this dataset
#' ?exampleData
#'
#' # Cluster data
#' # clustOutput <- mixGaussian::mixGaussianEM(dataset = exampleData,
#' #                                          membership = "none",
#' #                                          gmin = 1,
#' #                                          gmax = 5,
#' #                                          initMethod = "kmeans",
#' #                                          nInitIterations = 1)
#'
#' @author {Anjali Silva, \email{anjali@alumni.uoguelph.ca}}
#'
#' @references
#' Aitken, A. C. (1926). A series formula for the roots of algebraic and transcendental equations.
#' \emph{Proceedings of the Royal Society of Edinburgh}, 45, 14–22.
#'
#' B¨ohning, D., E. Dietz, R. Schaub, P. Schlattmann, and B. Lindsay (1994). The distribution
#' of the likelihood ratio for mixtures of densities from the one-parameter exponential family.
#' \emph{Annals of the Institute of Statistical Mathematics}, 46, 373–388.
#'
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
mixGaussianClassEM <- function(dataset,
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
    clusterResults[[gmodel]] <- mixGaussianClass(dataset = dataset,
                                                 initMethod = initMethod,
                                                 nInitIterations = nInitIterations,
                                                 membership = membership,
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

  class(RESULTS) <- "mixGaussianClass"
  return(RESULTS)

}

mixGaussianClass <- function(dataset,
                             G,
                             initMethod,
                             nInitIterations,
                             membership,
                             maxIterations) {

  dimensionality <- ncol(dataset)
  nObservations <- nrow(dataset)

  # If no intialization is requested by user
  if (nInitIterations == 0) {

    # basic initialization performed for parameters
    mu <- sigma <- num <- list()

    # kmeans initialization
    zValue <- matrix(0, ncol = G, nrow = nObservations)

    clsInd <- (membership == 0)
    for (i in 1:nObservations) {
      if(clsInd[i]) { # if no membership known
        zValue[i, ] <- 1 / G
      } else {
          zValue[i, membership[i]] <- 1
      }
    }
    z1 <- as.vector(t(zValue))

    piG <- colSums(zValue) / nObservations


    for (g in 1:G) {
        obs <- which(zValue[, g] == 1)
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



  # Start classification
  itOuter <- 1
  aloglik <- logLikelihood <- NULL
  conv <- aloglik[c(1, 2, 3)] <- 0

  while(! conv) {

    # Updating mu using known
    # membership values
    for(g in 1:G){
      mu[[g]] <- colSums(zValue[!clsInd, g] * dataset[!clsInd, ]) /
        sum(zValue[!clsInd, g])
    }

    # Updating covariance matrix
    # using known membership values
    for(g in 1:G){
      labelledObs <- c(1:nObservations)[!clsInd] # labelled observations
      for(i in labelledObs){ # only calculate for labelled observations
        num[[i]] <- zValue[i, g]*((dataset[i, ] - mu[[g]])%*%(t(dataset[i, ] - mu[[g]])))
      }
      sigma[[g]] <- Reduce('+', num)/sum(zValue[!clsInd, g])
    }



    # Update zvalue for unknown memberships
    # Temporary group membership (z) values
    tempZ <- matrix(0, length(which(clsInd == TRUE)), G)

    for(g in 1:G) {
      tempZ[which(clsInd == TRUE), g] <- piG[g] * mvtnorm::dmvnorm(x = dataset[clsInd, ],
                                              mean = mu[[g]],
                                              sigma = sigma[[g]])
    }

    # Calculate zValue value for unlabeled observations
    zValueUnobserved <- tempZ / rowSums(tempZ)
    zValue[which(clsInd == TRUE), ] <- zValueUnobserved


    # Calculate log-likelihood
    logLikelihood[itOuter] <- sum(log(rowSums(zValue)))


    # Update pig
    piG <- colSums(zValue) / nObservations

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

