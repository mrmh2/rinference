setwd('../')
source('nested.R')
source('models.R')

context("Testing prior generation functions")

test_that("generatePriorSamples works", {
  lower.bounds <- c(0, 3, 5)
  upper.bounds <- c(10, 6, 55)

  prior.samples <- generatePriorSamples(100, lower.bounds, upper.bounds)

  expect_equal(ncol(prior.samples), 100)
  expect_equal(nrow(prior.samples), 3)
  expect_true(all(prior.samples > lower.bounds))
  expect_true(all(prior.samples < lower.bounds))
})

context("Testing components")

test_that("makeStep works", {
  current.values <- c(3, 5, 7)
  step <- c(1, 10, 0.1)

  for (i in 1:10) {
    new.values <- makeStep(current.values, step)
    expect_false(any(new.values == current.values))
    within.bounds = abs(new.values - current.values) < step
    expect_true(all(within.bounds))
  }
})

## test_that("makeBoundedStep works", {
##   current.values <- c(3, 5, 7)
##   step <- c(5, 10, 20)
##   lower.bounds <- c(4, 6, 0)
##   upper.bounds <- c(7, 12, 20)

##   makeBoundedStep(current.values, step, lower.bounds, upper.bounds)
## })

test_that("makeBoundedSingleStep works", {
  current.value = 1
  step <- 0.5
  lower.bound <- 0.7
  upper.bound <- 1.3

  new.values <- replicate(100, makeBoundedSingleStep(current.value, step,
                                                     lower.bound, upper.bound))

  expect_true(all(new.values > lower.bound))
  expect_true(all(new.values < upper.bound))
})

test_that("makeMcmcStep works", {
  data <- read.csv('data/flowering-data.csv', sep=',', header=T)

  llh.sigma.factor <- 0.1
  linearModelLlFun <- function(params) { 
    return(generalLogLikelihood(linearModel, params, data, llh.sigma.factor))
  }

  current.values <- c(0.1, 0.1)

  lower.bounds <- c(0, 0)
  upper.bounds <- c(0.5, 0.5)

  steps <- 0.1 * (upper.bounds - lower.bounds)

  for (i in 1:10) {
    llMin <- linearModelLlFun(current.values)

    ret <- makeMcmcStep(current.values, linearModelLlFun, llMin, steps, lower.bounds, upper.bounds)

    if (any(ret$accepted) ) {
      expect_true(linearModelLlFun(ret$new.values) > llMin)
    }
    
    current.vales <- ret$new.values
  }

})  
  
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
  
  sigmoidal.params <- c(1, 10, 0.5, 0)
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
  expect_equal(logevidence, -109.1228, tolerance=0.01)

  ret <- nestedSampling(linearModelLlFun, prior.samples, bounds, 20)
  logevidence <- ret$logevidence
  expect_equal(logevidence, -206.944512, tolerance=0.01)
})
