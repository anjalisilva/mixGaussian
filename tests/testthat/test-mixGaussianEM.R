context("Checking for mixGaussian performance")
library(mixGaussian)

test_that("Checking clustering results", {

  # Generating simulated data
  G <- 2 # number of true clusters/components
  dimension <- 5
  nObservations <- 200
  piGTrue <- c(0.6, 0.4)

  set.seed(1234)
  mean1 <- rep(1, dimension)
  mean2 <- rep(6, dimension)
  sigma1 <- diag(dimension) * 0.5
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


  expect_that(length(clustOutput), equals(14))
  expect_that(clustOutput, is_a("mixGaussianEM"))
  expect_that(clustOutput$initalizationMethod, equals("kmeans"))
  numPara <- c(20, 41, 62)
  expect_that(clustOutput$numbParameters, equals(numPara))
  trueMembership <- c(rep(1, nObservations * piGTrue[1]),
                      rep(2, nObservations * piGTrue[2]))
  expect_that(clustOutput$trueLabels, equals(trueMembership))
  expect_that(trunc(clustOutput$ICLresults$ICLmodelselected), equals(2))
  expect_that(trunc(clustOutput$AICresults$AICmodelselected), equals(2))
  expect_that(trunc(clustOutput$AIC3results$AIC3modelselected), equals(2))
  expect_that(trunc(clustOutput$BICresults$BICmodelselected), equals(2))
})


context("Checking for invalid user input")
test_that("Data clustering error upon invalid user input", {

  # Generating simulated data
  G <- 2 # number of true clusters/components
  dimension <- 5
  nObservations <- 100
  piGTrue <- c(0.6, 0.4)

  set.seed(1234)
  mean1 <- rep(1, dimension)
  mean2 <- rep(6, dimension)
  sigma1 <- diag(dimension) * 0.5
  sigma2 <- diag(dimension) * 1

  component1 <- mvtnorm::rmvnorm(nObservations * piGTrue[1], mean = mean1, sigma = sigma1)
  component2 <- mvtnorm::rmvnorm(nObservations * piGTrue[2], mean = mean2, sigma = sigma2)
  dataset <- rbind(component1, component2)

  # Clustering
  clustOutput <- mixGaussian::mixGaussianEM(dataset = dataset,
                                            membership = c(rep(1, nObservations * piGTrue[1]),
                                                           rep(2, nObservations * piGTrue[2])),
                                            gmin = 1,
                                            gmax = 2,
                                            initMethod = "kmeans",
                                            nInitIterations = 1)

  # dataset provided as character
  expect_error(mixGaussianEM(dataset = "dataset",
                              membership = "none",
                              gmin = 1,
                              gmax = 2,
                              initMethod = "kmeans",
                              nInitIterations = 2))

  # dataset provided as logical
  expect_error(mixGaussianEM(dataset = NA,
                             membership = "none",
                             gmin = 1,
                             gmax = 2,
                             initMethod = "kmeans",
                             nInitIterations = 2))

  # Incorrect size for true membership
  trueMembership <- c(rep(1, nObservations * piGTrue[1]),
                      rep(2, nObservations * piGTrue[2]))
  expect_error(mixGaussianEM(dataset = dataset,
                             membership = trueMembership[-1],
                             gmin = 1,
                             gmax = 2,
                             initMethod = "kmeans",
                             nInitIterations = 2))


  # Incorrect class for true membership
  expect_error(mixGaussianEM(dataset = dataset,
                             membership = "trueMembership",
                             gmin = 1,
                             gmax = 2,
                             initMethod = "kmeans",
                             nInitIterations = 2))

  # Incorrect input type for gmin
  expect_error(mixGaussianEM(dataset = dataset,
                             membership = "none",
                             gmin = "1",
                             gmax = 2,
                             initMethod = "kmeans",
                             nInitIterations = 2))

  # Incorrect input for gmin and gmax
  expect_error(mixGaussianEM(dataset = dataset,
                             membership = "none",
                             gmin = 1,
                             gmax = "2",
                             initMethod = "kmeans",
                             nInitIterations = 2))



  # Incorrect input type for initMethod
  expect_error(mixGaussianEM(dataset = dataset,
                             membership = "none",
                             gmin = 1,
                             gmax = 2,
                             initMethod = "other",
                             nInitIterations = 2))

  # Incorrect input type for initMethod
  expect_error(mixGaussianEM(dataset = dataset,
                             membership = "none",
                             gmin = 1,
                             gmax = 2,
                             initMethod = NA,
                             nInitIterations = 2))

  # Incorrect input type for nInitIterations
  expect_error(mixGaussianEM(dataset = dataset,
                             membership = "none",
                             gmin = 1,
                             gmax = 2,
                             initMethod = "kmeans",
                             nInitIterations = "1"))

})
# [END]
