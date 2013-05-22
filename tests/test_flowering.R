setwd('../')
source('nested.R')
source('models.R')

context("Testing flowering time code")

test_that("linearModel works", {
  x <- seq(0, 1, 0.01)
  
  linear.params <- c(1, 0.5)
  y.linear <- linearModel(linear.params, x)
  y.exact <- seq(0.5, 1.5, 0.01)

  expect_equal(y.exact, y.linear)
})

test_that("sigmoidModel works", {
  x <- seq(0, 1, 0.01)
  
  sigmoidal.params <- c(1, 10, 0.5)
  y.sigmoidal <- sigmoidModel(sigmoidal.params, x)
  #plot(x, y.sigmoidal, type='l')

  expect_equal(y.sigmoidal[10], 0.0163025, tolerance=0.0001)
  expect_equal(y.sigmoidal[51], 0.5, tolerance=0.0001)
  expect_equal(y.sigmoidal[90], 0.9801597, tolerance=0.0001)
})

test_that("prior generation works", {
})

test_that("nestedSampling works", {

  set.seed(0)
  
  prior.size <- 10
  prior.lower.bounds <- c(0, 0)
  prior.upper.bounds <- c(0.5, 0.001)
  prior.samples <- replicate(prior.size, generateAPrior(prior.lower.bounds, prior.upper.bounds))

  data <- read.csv('data/flowering-data.csv', sep=',', header=T)

  llh.sigma.factor <- 0.1
  linearModelLlFun <- function(params) { 
    return(generalLogLikelihood(linearModel, params, data, llh.sigma.factor))
  }

  lower.bounds <- prior.lower.bounds
  upper.bounds <- prior.upper.bounds
  bounds <- cbind(lower.bounds, upper.bounds)
  
  steps <- c(0.01, 0.0001)
  ret <- nestedSampling(linearModelLlFun, prior.samples, bounds, 20, steps=steps)
  logevidence <- ret$logevidence
  expect_equal(logevidence, -134.0461, tolerance=0.01)

  ret <- nestedSampling(linearModelLlFun, prior.samples, bounds, 20)
  logevidence <- ret$logevidence
  expect_equal(logevidence, -101.6245, tolerance=0.01)
})
