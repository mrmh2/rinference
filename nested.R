########################################################################
# nested.R - code to estimate posterior and evidence via nested sampling
# (c) Matthew Hartley, Nick Pullen and Richard Morris
########################################################################

logPlus <- function(a, b) {
  # Sum a and b via their logarithms, such that if:
  # a = log(x) and b = log(y)
  # logPlus(a, b) returns log(x + y)
  if (a > b) {
      sum = a + log(1+exp(b-a))
  } else {
      sum = b + log(1+exp(a-b))
  }
  return(sum)
}

makeStep <- function(current.values, step) {
  # Make random step in parameter space
  return(current.values + step * runif(1, -1, 1))
}

makeBoundedStep <- function(current.values, step, lower.bounds, upper.bounds) {
  # Make step in parameter space. If step exceeds bounds, instead pick point in
  # space randomly chosen using uniform distribution.

  step.attempt <- makeStep(current.values, step)
  # Generate vector of booleans that are true where we exceed bounds
  outside.bounds <- (step.attempt < lower.bounds) + (step.attempt > upper.bounds) > 0
  # Replace values outside the bounds by randomly chosen point in parameter space
  step.attempt[which(outside.bounds)] = runif(1, lower.bounds, upper.bounds)

  return(step.attempt)
}

makeMcmcStep <- function(current.values, llFun, llMin, steps, lower.bounds, upper.bounds) {
  # Make single MCMC step from current.values as follows:
  # 1) Attempt step in each parameter direction
  # 2) If step is successful (ll bigger than llMin), note this
  # 
  # Returns:
  #   A list containing:
  #    accepted: a vector of boolean values describing whether each parameter's
  #              step was accepted
  #    new.values: a vector of the new parameter values after the step

  n <- length(current.values)
  accepted <- rep(FALSE, n)
  for (i in 1:n) {
    empty <- rep(0, n)
    empty[i] <- steps[i]
    # Attempt step along one axis in parameter space
    candidate <- makeBoundedStep(current.values, empty, lower.bounds, upper.bounds)
    if(llFun(candidate) > llMin) {
      current.values <- candidate
      accepted[i] = TRUE
    }
  }
  return(list(accepted=accepted, new.values=current.values))
}

explore <- function(current.values, steps, llMin, llFun, lower.bounds, upper.bounds) {
  # Explore parameter space around supplied point
  # returns new point and new step size
  msub <- 20

  accepted <- vector()
  for(i in 1:msub) {
    ret <- makeMcmcStep(current.values, llFun, llMin, steps, lower.bounds, upper.bounds)
    accepted <- cbind(accepted, ret$accepted)
    current.values <- ret$new.values
  }

  # DEBUG
  #cat(sum(accepted[1,]), " ", sum(accepted[2,]), "\n")

  for(i in 1:length(current.values)) {
    ratio <- sum(accepted[i,]) / length(accepted[i,])
    if (ratio > 0.6) {
        steps[i] <- steps[i] * (1 + 2 * (ratio - 0.6) / 0.4)
    }
    if (ratio < 0.4) {
        steps[i] <- steps[i] / (1 + 2 * ((0.4 - ratio) / 0.4))
    }
    steps[i] <- min(steps[i], 0.1 * (upper.bounds[i] - lower.bounds[i]))

    #cat(sprintf("%f, %f\n", ratio, step[i]))
  }

  return(list(new.values=current.values, new.step=steps, new.ll=llFun(current.values)))
}

prior <- function(  xmin, xmax, n) {
# Return a vector of length n samples from the uniform distribution in the
# range (xmin, xmax)

    return( runif(n, xmin, xmax) )
}

generatePriorSamples <- function(n.samples, bounds)  {

  prior.lower.bounds <- bounds[,1]
  prior.upper.bounds <- bounds[,2]

  prior.samples <- replicate(n.samples,  generateAPrior(prior.lower.bounds, prior.upper.bounds))

  return(prior.samples)
}

generateAPrior <- function(prior.lower.bounds, prior.upper.bounds) {
  # Generate a single prior drawn from the uniform distribution
  # 
  # Args:
  #   prior.lower.bounds - a vector of lower bounds
  #   prior.upper.bounds - a vector of upper bounds
  # 
  # Returns:
  #   A vector of values drawn from the uniform distribution with the
  #   bounds specified

  n <- length(prior.upper.bounds)
  prior <- numeric(n)

  for (i in 1:n) {
    prior[i] <- runif(1, prior.lower.bounds[i], prior.upper.bounds[i])
  }

  return(prior)
}

calculateEvidence <- function(posterior.samples, prior.size) {
  # Calculates the evidence from the set of posterior samples with
  # loglikelihoods.
  #
  # Args:
  #  posterior.samples: The posterior samples, together with log likelihoods 
  #                     as an m x n matrix, where m is the number of parameters
  #                     plus one and n is the number of posterior samples. The
  #                     first row of the matrix should be the log likelihoods.
  #  prior.size:        The numbers of samples in the prior.
  #
  # Returns:
  #   The evidence, a scalar.

  # FIXME - problems with machine precision for some reason?
  #logZ <- -.Machine$double.xmin
  posterior.size <- ncol(posterior.samples)
  logZ <- -2.2e308
  Xlast <- 1
  for (i in 1:(posterior.size-prior.size)) {
    lL <- posterior.samples[1,i]
    X <- exp(-i / prior.size)
    lw <- log(Xlast - X)

    logZ <- logPlus(logZ, lw + lL)

    Xlast <- X    
  }

  for (i in 1:prior.size) {
    lL <- posterior.samples[1, i + posterior.size - prior.size]
    logZ <- logPlus(logZ, lw + lL)
  }

  return(logZ)
}

calculatePosterior <- function(posterior.size, ordered.samples, steps, llFun,
                               lower.bounds, upper.bounds) {
  # Generate the posterior samples for the given model.

  # Initialise an empty matrix to hold posterior samples & their loglikelihoods
  prior.size = ncol(ordered.samples)
  posterior.samples <- matrix(nrow=nrow(ordered.samples))
  first <- TRUE

  while ( ncol(posterior.samples) < posterior.size ) {
    # Add the worst sample to the posterior samples
    if ( first ) {
      posterior.samples[,1] <- ordered.samples[,1]
      first <- FALSE
    } else {
      posterior.samples <- cbind(posterior.samples, ordered.samples[,1])
    }

    # Randomly select a sample that isn't the worst
    selected.point <- sample(2:prior.size, 1)

    # Explore around that point, updating step size as we go
    llMin <- min(ordered.samples[1,])
    ret <- explore(ordered.samples[2:nrow(ordered.samples),selected.point], 
           steps, llMin, llFun, lower.bounds, upper.bounds)
    steps <- ret$new.step
    new.point <- ret$new.values
    new.ll <- ret$new.ll

    # Replace the worst sample by the results of explore()
    ordered.samples[,1] <- c(new.ll, new.point)

    # Re-sort the samples (not very efficient, but we don't care)
    ordered.samples <- ordered.samples[,order(ordered.samples[1,])]
  }

  return(posterior.samples)
}

nestedSampling <- function(llFun, prior.samples, bounds, posterior.size, steps=NULL) {
  # Perform nested sampling for the given model (expressed through the log
  # likelihood function).
  #
  # Args:
  #   llFun:          The loglikelihood function, which should be of the form:
  #                   llFun(params) where params is a vector containing a single
  #                   instance of parameters, and should return the log 
  #                   likelihood for those parameters.
  #   prior.samples:  The samples from the prior, as a m x n matrix, where m is
  #                   the number of parameters + 1 and n is the size of the 
  #                   prior.
  #   bounds:         The bounds of the parameter values, as a m x 2 matrix, 
  #                   where m is the number of parameters, the first column is 
  #                   the lower bounds and .the second the upper bounds
  #   posterior.size: The number of desired sample points in the posterior to
  #                   be calculated.
  #   steps:          A vector of the initial step sizes to use for exploring
  #                   the parameter space, of size m where m is the number of
  #                   parameters.
  #
  # Returns:
  #   logevidence:       The logarithm of the evidence - P(model|data).
  #   posterior.samples: The samples from the posterior, together with their
  #                      log likelihoods as a m x n matrix, where m is the
  #                      number of parameters + 1 and n is the number of
  #                      posterior samples. The log likelhood values are the
  #                      top row of this matrix.

  prior.size = ncol(prior.samples)
  
  # Calculate the log likelihoods for the prior samples
  ll.values <- apply(prior.samples, 2, llFun)
  evaluated.samples = rbind(ll.values, prior.samples)
  ordered.samples <- evaluated.samples[,order(evaluated.samples[1,])]

  lower.bounds <- bounds[,1]
  upper.bounds <- bounds[,2]

  # If steps is not set calculate from range of bounds
  if (is.null(steps)) {
    steps <- (upper.bounds - lower.bounds) / 100
  }

  # Generate the posterior samples
  posterior.samples <- calculatePosterior(posterior.size, ordered.samples,
    steps, llFun, lower.bounds, upper.bounds)

  # Calculate the evidence from the posterior samples
  logZ <- calculateEvidence(posterior.samples, prior.size)

  return(list(logevidence=unname(logZ), posterior=posterior.samples))

}
