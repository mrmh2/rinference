setwd('../')
source('flowering.R')

context("Testing flowering time code")

test_that("calculateEvidence works", {

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
  
  ret <- calculateEvidence(linearModelLlFun, 20, prior.samples)

  logevidence <- ret$logevidence

  expect_equal(logevidence, -134.0461, tolerance=0.01)
})
