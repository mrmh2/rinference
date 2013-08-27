library(ggplot2)
# Load the nested sampling code.
source('nested.R')
# Load the models
source('models.R')

#############################################################
# Example 1 - fitting a linear model to flowering time data #
#############################################################

## # Load the flowering time data
## #data <- read.csv('data/flowering-data.csv', sep=',', header=T)
## data <- read.table('data/TFL1-FT-data', header=FALSE)

## # Generate the log likelihood function for a linear model for our data
## llh.sigma.factor <- 0.1
## linearModelLlFun <- function(params) { 
##   return(generalLogLikelihood(linearModel, params, data, llh.sigma.factor))
## }

## gradLinearModel <- function(params) {
##   return(gradLogLikelihood(linearModel, params, data, llh.sigma.factor))
## }

## # Set the prior size and bounds, then generate a set of prior samples
## prior.size <- 25
## #prior.bounds <- matrix(c(0, 0, 0.2, 0.025), nrow=2, ncol=2) #test bounds to show problem, don't use otherwise!
## prior.bounds <- matrix(c(0, 0, 0.5, 0.5), nrow=2, ncol=2)
## #prior.bounds <- matrix(c(0, -2., 10, 10), nrow=2, ncol=2) 
## prior.samples <- generatePriorSamples(prior.size, prior.bounds)

## # Set the bounds for the parameter exploration
## bounds <- prior.bounds

## # Set the number of posterior samples we'd like to obtain
## posterior.samples <- 200

## # Call the nested sampling routine - this returns the posterior and log evidence
## #ret <- nestedSampling(linearModelLlFun, prior.samples, bounds, posterior.samples)
## #ret <- nestedSampling(linearModelLlFun, prior.samples, bounds, posterior.samples, mcmcMethod='HMC')
## #ret <- nestedSampling(linearModelLlFun, prior.samples, bounds, posterior.samples, gradLinearModel)
## ret <- nestedSampling(linearModelLlFun, prior.samples, bounds, posterior.samples, gradLinearModel, mcmcMethod='HMC')
## posterior <- ret$posterior

## # Plot the data, along with the linear model generated 
## #plot(data)
## #abline(posterior[3, posterior.samples], posterior[2, posterior.samples])

## # Display the evidence
## cat('log evidence = ', ret$logevidence, '\n')

## # Find the best set of parameters
## best.params <- posterior[2:3, posterior.samples]
## cat('best NS params =', rev(best.params), '\n')

## # Find the exact parameters via R's linear model fitting function
## lm.results <- lm(as.matrix(data[2]) ~ as.matrix(data[1]))
## lm.coefficients <- unname(lm.results$coefficients)
## cat('lm coefficients=',lm.coefficients, '\n')

## # Calculate the percentage error in our parameter estimates
## param.error <-(rev(best.params) - lm.coefficients) / lm.coefficients
## cat('parameter error vector = ', param.error, '\n')

  #################################################################
  ## Example 2 - fitting a sigmoidal model to flowering time data #
  #################################################################
#set.seed(1)
  #
  ## Load the flowering time data
  #data <- read.csv('data/flowering-data.csv', sep=',', header=T)
  data <- read.table('data/TFL1-FT-data', header=FALSE)

  ## Generate the log likelihood function for a sigmoidal model for our data
  llh.sigma.factor <- 0.1
  sigmoidalModelLlFun <- function(params) { 
    return(generalLogLikelihood(sigmoidModel, params, data, llh.sigma.factor))
  }

  gradSigmoidalModel <- function(params) {
    return(gradLogLikelihood(sigmoidModel, params, data, llh.sigma.factor))
  }
  #
  ## Set the prior size and bounds, then generate a set of prior samples
  prior.size <- 20
  #prior.lower.bounds <- c(0, 0, 1, 0)
  #prior.upper.bounds <- c(5., 5., 9., 5.)
  prior.lower.bounds <- c(0.5, 0, 2, 0)
  prior.upper.bounds <- c(1.5, 4., 8., 0.5)
  prior.bounds = cbind(prior.lower.bounds, prior.upper.bounds)
  prior.samples <- generatePriorSamples(prior.size, prior.bounds)
  #
  ## Set the bounds for the parameter exploration
  bounds <- prior.bounds
  #
  ## Set the number of posterior samples we'd like to obtain
  posterior.samples <- 200
  #
  ## Call the nested sampling routine - this returns the posterior and log evidence
  ret <- nestedSampling(sigmoidalModelLlFun, prior.samples, bounds, posterior.samples, gradient=gradSigmoidalModel, mcmcMethod='HMC')
#  ret <- nestedSampling(sigmoidalModelLlFun, prior.samples, bounds, posterior.samples, gradient=gradSigmoidalModel)
  posterior <- ret$posterior
  #
  ## Plot the data, along with the sigmoidal fit
  x <- as.matrix(data[1])
  best.params <- posterior[2:5, posterior.samples]
  cat('best NS params =', best.params, '\n')
  y.sigmoidal <- sigmoidModel(best.params, x)
  y.actual <- as.matrix(data[2])
  plot(x, y.actual)
  lines(x, y.sigmoidal,lty=1)
  #
  ## Display the evidence
  cat('log evidence = ', ret$logevidence, '\n')
  #
  ## Calculate fit using R's fitting function
  #
  #sigfit <- nls(y.actual ~ k4 + (k1 - k4) / (1 + exp(-k2 * (x - k3))),
  #            start=list(k1=1.4, k2=2, k3=5.5, k4=0.4))
  #
  #fit.params <- coef(sigfit)
  #print(fit.params)
  #y.fit <- sigmoidModel(fit.params, x)
  #lines(x, y.fit,lty=2)
#
#
#
