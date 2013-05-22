library(ggplot2)
# Load the nested sampling code.
source('nested.R')
# Load the models
source('models.R')

#############################################################
# Example 1 - fitting a linear model to flowering time data #
#############################################################

# Load the flowering time data
data <- read.csv('data/flowering-data.csv', sep=',', header=T)

# Generate the log likelihood function for a linear model for our data
llh.sigma.factor <- 0.1
linearModelLlFun <- function(params) { 
  return(generalLogLikelihood(linearModel, params, data, llh.sigma.factor))
}

# Set the prior size and bounds, then generate a set of prior samples
prior.size <- 20
prior.bounds <- matrix(c(0, 0, 0.5, 0.001), nrow=2, ncol=2)
prior.samples <- generatePriorSamples(prior.size, prior.bounds)

# Set the bounds for the parameter exploration
bounds <- prior.bounds

# Set the number of posterior samples we'd like to obtain
posterior.samples <- 200

# Call the nested sampling routine - this returns the posterior and log evidence
ret <- nestedSampling(linearModelLlFun, prior.samples, bounds, posterior.samples)
posterior <- ret$posterior

# Plot the data, along with the linear model generated 
plot(data)
abline(posterior[3, posterior.samples], posterior[2, posterior.samples])

# Display the evidence
cat('log evidence = ', ret$logevidence, '\n')

# Find the best set of parameters
best.params <- posterior[2:3, posterior.samples]

# Find the exact parameters via R's linear model fitting function
lm.results <- lm(as.matrix(data[2]) ~ as.matrix(data[1]))
lm.coefficients <- unname(lm.results$coefficients)

# Calculate the percentage error in our parameter estimates
param.error <-(rev(best.params) - lm.coefficients) / lm.coefficients
cat('parameter error vector = ', param.error, '\n')

################################################################
# Example 2 - fitting a sigmoidal model to flowering time data #
################################################################

# Load the flowering time data
data <- read.csv('data/flowering-data.csv', sep=',', header=T)

# Generate the log likelihood function for a sigmoidal model for our data
llh.sigma.factor <- 0.1
sigmoidalModelLlFun <- function(params) { 
  return(generalLogLikelihood(sigmoidModel, params, data, llh.sigma.factor))
}

# Set the prior size and bounds, then generate a set of prior samples
prior.size <- 20
prior.lower.bounds <- c(0.001, 500, 0.005)
prior.upper.bounds <- c(0.002, 1500, 0.01)
prior.bounds = cbind(prior.lower.bounds, prior.upper.bounds)
prior.samples <- generatePriorSamples(prior.size, prior.bounds)

# Set the bounds for the parameter exploration
bounds <- prior.bounds

# Set the number of posterior samples we'd like to obtain
posterior.samples <- 200

# Call the nested sampling routine - this returns the posterior and log evidence
ret <- nestedSampling(sigmoidalModelLlFun, prior.samples, bounds, posterior.samples)
posterior <- ret$posterior

# Plot the data, along with the sigmoidal fit
x <- as.matrix(data[1])
best.params <- posterior[2:4, posterior.samples]
y.sigmoidal <- sigmoidModel(best.params, x)
y.actual <- as.matrix(data[2])
plot(x, y.actual)
par(new=T)
plot(x, y.sigmoidal, type='l')

# Display the evidence
cat('log evidence = ', ret$logevidence, '\n')


