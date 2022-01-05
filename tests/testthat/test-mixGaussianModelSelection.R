# context("Checking for model selection performance")
# library(mixGaussian)

test_that("Checking AIC model selection", {

 # Generating simulated data
 G <- 2 # number of true clusters/components
 dimension <- 6
 nObservations <- 200
 piGTrue <- c(0.7, 0.3)

 set.seed(1234)
 mean1 <- rep(1, dimension)
 mean2 <- rep(4, dimension)
 sigma1 <- diag(dimension) * 1
 sigma2 <- diag(dimension) * 1

 component1 <- mvtnorm::rmvnorm(nObservations * piGTrue[1], mean = mean1, sigma = sigma1)
 component2 <- mvtnorm::rmvnorm(nObservations * piGTrue[2], mean = mean2, sigma = sigma2)
 dataset <- rbind(component1, component2)

 # Clustering
 clustOutput <- mixGaussian::mixGaussianEM(dataset = dataset,
                                            membership = c(rep(1, nObservations * piGTrue[1]),
                                                           rep(2, nObservations * piGTrue[2])),
                                            gmin = 1,
                                            gmax = 3,
                                            initMethod = "kmeans",
                                            nInitIterations = 1)
 # Model Selection
 AICmodel <- mixGaussian::AICFunction(logLikelihood = clustOutput$logLikelihood,
                                      nParameters = clustOutput$numbParameters,
                                      clusterRunOutput = clustOutput$allResults,
                                      gmin = clustOutput$gmin,
                                      gmax = clustOutput$gmax)

  expect_that(length(AICmodel), equals(4))
  expect_that(AICmodel, is_a("AIC"))
  expect_that(length(unique(AICmodel$AICmodelSelectedLabels)), equals(2))
  expect_that(AICmodel$allAICvalues, is_a("numeric"))
  expect_that(trunc(AICmodel$AICmodelselected), equals(2))
 })
context("Checking for invalid user input for AIC")
test_that("AIC model selection error upon invalid user input", {

  # Generating simulated data
  G <- 2 # number of true clusters/components
  dimension <- 6
  nObservations <- 200
  piGTrue <- c(0.7, 0.3)

  set.seed(1234)
  mean1 <- rep(1, dimension)
  mean2 <- rep(4, dimension)
  sigma1 <- diag(dimension) * 1
  sigma2 <- diag(dimension) * 1

  component1 <- mvtnorm::rmvnorm(nObservations * piGTrue[1], mean = mean1, sigma = sigma1)
  component2 <- mvtnorm::rmvnorm(nObservations * piGTrue[2], mean = mean2, sigma = sigma2)
  dataset <- rbind(component1, component2)

  # Clustering
  clustOutput <- mixGaussian::mixGaussianEM(dataset = dataset,
                                            membership = c(rep(1, nObservations * piGTrue[1]),
                                                           rep(2, nObservations * piGTrue[2])),
                                            gmin = 1,
                                            gmax = 3,
                                            initMethod = "kmeans",
                                            nInitIterations = 1)

  # logLikelihood provided as character
  expect_error(AICFunction(logLikelihood = "clustOutput$logLikelihood",
                          nParameters = clustOutput$numbParameters,
                          clusterRunOutput = clustOutput$allResults,
                          gmin = clustOutput$gmin,
                          gmax = clustOutput$gmax))

  # nParameters provided as character
  expect_error(AICFunction(logLikelihood = clustOutput$logLikelihood,
                           nParameters = "clustOutput$numbParameters",
                           clusterRunOutput = clustOutput$allResults,
                           gmin = clustOutput$gmin,
                           gmax = clustOutput$gmax))

  # gmin provided as character
  expect_error(AICFunction(logLikelihood = clustOutput$logLikelihood,
                          nParameters = clustOutput$numbParameters,
                          clusterRunOutput = clustOutput$allResults,
                          gmin = "clustOutput$gmin",
                          gmax = clustOutput$gmax))

  # gmax provided as character
  expect_error(AICFunction(logLikelihood = clustOutput$logLikelihood,
                            nParameters = clustOutput$numbParameters,
                            clusterRunOutput = clustOutput$allResults,
                            gmin = clustOutput$gmin,
                            gmax = "clustOutput$gmax"))

})


test_that("Checking AIC3 model selection", {

 # Generating simulated data
 G <- 2 # number of true clusters/components
 dimension <- 6
 nObservations <- 200
 piGTrue <- c(0.7, 0.3)

 set.seed(1234)
 mean1 <- rep(1, dimension)
 mean2 <- rep(4, dimension)
 sigma1 <- diag(dimension) * 1
 sigma2 <- diag(dimension) * 1

 component1 <- mvtnorm::rmvnorm(nObservations * piGTrue[1], mean = mean1, sigma = sigma1)
 component2 <- mvtnorm::rmvnorm(nObservations * piGTrue[2], mean = mean2, sigma = sigma2)
 dataset <- rbind(component1, component2)

 # Clustering
 clustOutput <- mixGaussian::mixGaussianEM(dataset = dataset,
                                            membership = c(rep(1, nObservations * piGTrue[1]),
                                                           rep(2, nObservations * piGTrue[2])),
                                            gmin = 1,
                                            gmax = 3,
                                            initMethod = "kmeans",
                                            nInitIterations = 1)
 # Model Selection
 AIC3model <- AIC3Function(logLikelihood = clustOutput$logLikelihood,
                           nParameters = clustOutput$numbParameters,
                           clusterRunOutput = clustOutput$allResults,
                           gmin = clustOutput$gmin,
                           gmax = clustOutput$gmax)

 expect_that(length(AIC3model), equals(4))
 expect_that(AIC3model, is_a("AIC3"))
 expect_that(length(unique(AIC3model$AIC3modelSelectedLabels)), equals(2))
 expect_that(AIC3model$allAIC3values, is_a("numeric"))
 expect_that(trunc(AIC3model$AIC3modelselected), equals(2))
})
context("Checking for invalid user input for AIC3")
test_that("AIC3 model selection error upon invalid user input", {

  # Generating simulated data
  G <- 2 # number of true clusters/components
  dimension <- 6
  nObservations <- 200
  piGTrue <- c(0.7, 0.3)

  set.seed(1234)
  mean1 <- rep(1, dimension)
  mean2 <- rep(4, dimension)
  sigma1 <- diag(dimension) * 1
  sigma2 <- diag(dimension) * 1

  component1 <- mvtnorm::rmvnorm(nObservations * piGTrue[1], mean = mean1, sigma = sigma1)
  component2 <- mvtnorm::rmvnorm(nObservations * piGTrue[2], mean = mean2, sigma = sigma2)
  dataset <- rbind(component1, component2)

  # Clustering
  clustOutput <- mixGaussian::mixGaussianEM(dataset = dataset,
                                            membership = c(rep(1, nObservations * piGTrue[1]),
                                                           rep(2, nObservations * piGTrue[2])),
                                            gmin = 1,
                                            gmax = 3,
                                            initMethod = "kmeans",
                                            nInitIterations = 1)

  # logLikelihood provided as character
  expect_error(AIC3Function(logLikelihood = "clustOutput$logLikelihood",
                            nParameters = clustOutput$numbParameters,
                            clusterRunOutput = clustOutput$allResults,
                            gmin = clustOutput$gmin,
                            gmax = clustOutput$gmax))

  # nParameters provided as character
  expect_error(AIC3Function(logLikelihood = clustOutput$logLikelihood,
                            nParameters = "clustOutput$numbParameters",
                            clusterRunOutput = clustOutput$allResults,
                            gmin = clustOutput$gmin,
                            gmax = clustOutput$gmax))

  # gmin provided as character
  expect_error(AIC3Function(logLikelihood = clustOutput$logLikelihood,
                            nParameters = clustOutput$numbParameters,
                            clusterRunOutput = clustOutput$allResults,
                            gmin = "clustOutput$gmin",
                            gmax = clustOutput$gmax))

  # gmax provided as character
  expect_error(AIC3Function(logLikelihood = clustOutput$logLikelihood,
                            nParameters = clustOutput$numbParameters,
                            clusterRunOutput = clustOutput$allResults,
                            gmin = clustOutput$gmin,
                            gmax = "clustOutput$gmax"))

})


test_that("Checking BIC model selection", {

  # Generating simulated data
  G <- 2 # number of true clusters/components
  dimension <- 6
  nObservations <- 200
  piGTrue <- c(0.7, 0.3)

  set.seed(1234)
  mean1 <- rep(1, dimension)
  mean2 <- rep(4, dimension)
  sigma1 <- diag(dimension) * 1
  sigma2 <- diag(dimension) * 1

  component1 <- mvtnorm::rmvnorm(nObservations * piGTrue[1], mean = mean1, sigma = sigma1)
  component2 <- mvtnorm::rmvnorm(nObservations * piGTrue[2], mean = mean2, sigma = sigma2)
  dataset <- rbind(component1, component2)

  # Clustering
  clustOutput <- mixGaussian::mixGaussianEM(dataset = dataset,
                                            membership = c(rep(1, nObservations * piGTrue[1]),
                                                           rep(2, nObservations * piGTrue[2])),
                                            gmin = 1,
                                            gmax = 3,
                                            initMethod = "kmeans",
                                            nInitIterations = 1)

  BICmodel <- BICFunction(logLikelihood = clustOutput$logLikelihood,
                          nParameters = clustOutput$numbParameters,
                          nObservations = nrow(clustOutput$dataset),
                          clusterRunOutput = clustOutput$allResults,
                          gmin = clustOutput$gmin,
                          gmax = clustOutput$gmax)

  expect_that(length(BICmodel), equals(4))
  expect_that(BICmodel, is_a("BIC"))
  expect_that(length(unique(BICmodel$BICmodelSelectedLabels)), equals(2))
  expect_that(BICmodel$allBICvalues, is_a("numeric"))
  expect_that(trunc(BICmodel$BICmodelselected), equals(2))
})

context("Checking for invalid user input for BIC")
test_that("BIC model selection error upon invalid user input", {


  # logLikelihood provided as character
  expect_error(BICFunction(logLikelihood = "clustOutput$logLikelihood",
                           nParameters = clustOutput$numbParameters,
                           nObservations = nrow(clustOutput$dataset),
                           clusterRunOutput = clustOutput$allResults,
                           gmin = clustOutput$gmin,
                           gmax = clustOutput$gmax))

  # nParameters provided as character
  expect_error(BICFunction(logLikelihood = clustOutput$logLikelihood,
                          nParameters = "clustOutput$numbParameters",
                          nObservations = nrow(clustOutput$dataset),
                          clusterRunOutput = clustOutput$allResults,
                          gmin = clustOutput$gmin,
                          gmax = clustOutput$gmax))

  # gmin provided as character
  expect_error(BICFunction(logLikelihood = clustOutput$logLikelihood,
                          nParameters = clustOutput$numbParameters,
                          nObservations = nrow(clustOutput$dataset),
                          clusterRunOutput = clustOutput$allResults,
                          gmin = "clustOutput$gmin",
                          gmax = clustOutput$gmax))

  # gmax provided as character
  expect_error(BICFunction(logLikelihood = clustOutput$logLikelihood,
                          nParameters = clustOutput$numbParameters,
                          nObservations = nrow(clustOutput$dataset),
                          clusterRunOutput = clustOutput$allResults,
                          gmin = clustOutput$gmin,
                          gmax = "clustOutput$gmax"))
})


test_that("Checking ICL model selection", {

  # Generating simulated data
  G <- 2 # number of true clusters/components
  dimension <- 6
  nObservations <- 200
  piGTrue <- c(0.7, 0.3)

  set.seed(1234)
  mean1 <- rep(1, dimension)
  mean2 <- rep(4, dimension)
  sigma1 <- diag(dimension) * 1
  sigma2 <- diag(dimension) * 1

  component1 <- mvtnorm::rmvnorm(nObservations * piGTrue[1], mean = mean1, sigma = sigma1)
  component2 <- mvtnorm::rmvnorm(nObservations * piGTrue[2], mean = mean2, sigma = sigma2)
  dataset <- rbind(component1, component2)

  # Clustering
  clustOutput <- mixGaussian::mixGaussianEM(dataset = dataset,
                                            membership = c(rep(1, nObservations * piGTrue[1]),
                                                           rep(2, nObservations * piGTrue[2])),
                                            gmin = 1,
                                            gmax = 3,
                                            initMethod = "kmeans",
                                            nInitIterations = 1)

  ICLmodel <- ICLFunction(logLikelihood = clustOutput$logLikelihood,
                          nParameters = clustOutput$numbParameters,
                          nObservations = nrow(clustOutput$dataset),
                          clusterRunOutput = clustOutput$allResults,
                          gmin = clustOutput$gmin,
                          gmax = clustOutput$gmax)

 expect_that(length(ICLmodel), equals(4))
 expect_that(ICLmodel, is_a("ICL"))
 expect_that(length(unique(ICLmodel$ICLmodelSelectedLabels)), equals(2))
 expect_that(ICLmodel$allICLvalues, is_a("numeric"))
 expect_that(trunc(ICLmodel$ICLmodelselected), equals(2))
})


context("Checking for invalid user input for ICL")
test_that("ICL model selection error upon invalid user input", {

  # Generating simulated data
  G <- 2 # number of true clusters/components
  dimension <- 6
  nObservations <- 200
  piGTrue <- c(0.7, 0.3)

  set.seed(1234)
  mean1 <- rep(1, dimension)
  mean2 <- rep(4, dimension)
  sigma1 <- diag(dimension) * 1
  sigma2 <- diag(dimension) * 1

  component1 <- mvtnorm::rmvnorm(nObservations * piGTrue[1], mean = mean1, sigma = sigma1)
  component2 <- mvtnorm::rmvnorm(nObservations * piGTrue[2], mean = mean2, sigma = sigma2)
  dataset <- rbind(component1, component2)

  # Clustering
  clustOutput <- mixGaussian::mixGaussianEM(dataset = dataset,
                                            membership = c(rep(1, nObservations * piGTrue[1]),
                                                           rep(2, nObservations * piGTrue[2])),
                                            gmin = 1,
                                            gmax = 3,
                                            initMethod = "kmeans",
                                            nInitIterations = 1)

  expect_error(ICLFunction(logLikelihood = "clustOutput$logLikelihood",
                            nParameters = clustOutput$numbParameters,
                            nObservations = nrow(clustOutput$dataset),
                            clusterRunOutput = clustOutput$allResults,
                            gmin = clustOutput$gmin,
                            gmax = clustOutput$gmax))

  # nParameters provided as character
  expect_error(ICLFunction(logLikelihood = mplnResults$logLikelihood,
                          nParameters = "mplnResults$numbParameters",
                          nObservations = nrow(mplnResults$dataset),
                          clusterRunOutput = mplnResults$allResults,
                          gmin = mplnResults$gmin,
                          gmax = mplnResults$gmax,
                          parallel = FALSE))

  # gmin provided as character
  expect_error(ICLFunction(logLikelihood = clustOutput$logLikelihood,
                            nParameters = clustOutput$numbParameters,
                            nObservations = nrow(clustOutput$dataset),
                            clusterRunOutput = clustOutput$allResults,
                            gmin = "clustOutput$gmin",
                            gmax = clustOutput$gmax))

  # gmax provided as character
  expect_error(ICLFunction(logLikelihood = clustOutput$logLikelihood,
                            nParameters = clustOutput$numbParameters,
                            nObservations = nrow(clustOutput$dataset),
                            clusterRunOutput = clustOutput$allResults,
                            gmin = clustOutput$gmin,
                            gmax = "clustOutput$gmax"))

})
# [END]
