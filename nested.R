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

calculateEvidence <- function(llFun, posterior.size, prior.samples) {
  # Args:
  #  llFun - the loglikelihood function, which should be of the form:
  #          llFun(params) where params is a vector containing a single
  #          instance of parameters, and should return the log likelihood
  #          for those parameters
  #

  prior.size = dim(prior.samples)[2]
  
  ll.values <- apply(prior.samples, 2, llFun)
  evaluated.samples = rbind(ll.values, prior.samples)
  ordered.samples <- evaluated.samples[,order(evaluated.samples[1,])]

  lower.bounds <- c(0, 0)
  upper.bounds <- c(0.5, 0.001)
  steps <- c(0.01, 0.0001)

  # Initialise an empty matrix to hold posterior samples & their loglikelihoods
  posterior.samples <- matrix(nrow=3)
  first <- TRUE

  while ( dim(posterior.samples)[2] < posterior.size ) {
    # Add the worst sample to the posterior samples
    if ( first ) {
      posterior.samples[,1] <- ordered.samples[,1]
      first <- FALSE
    } else {
      posterior.samples <- cbind(posterior.samples, ordered.samples[,1])
    }

    # Randomly select a sample that isn't the worst
    selected.point <- sample(2:prior.size, 1)

    # Explore around that point
    llMin <- min(ordered.samples[1,])
    ret <- explore(ordered.samples[2:3,selected.point], steps, llMin, llFun, lower.bounds, upper.bounds)
    steps <- ret$new.step
    new.point <- ret$new.values
    new.ll <- ret$new.ll

    # Replace the worst sample by the results of explore()
    ordered.samples[,1] <- c(new.ll, new.point)

    # Re-sort the samples (not very efficient, but we don't care)
    ordered.samples <- ordered.samples[,order(ordered.samples[1,])]
  }

  #logZ <- .Machine$double.xmin
  logZ <- -1e300
  Xlast <- 1
  for (i in 1:posterior.size) {
    lL <- posterior.samples[1,i]
    X <- exp(-i / prior.size)
    lw <- log(Xlast - X)
    #cat(lw, ' ', lL, ' ', logZ + lw + lL, '\n')

    logZ <- logPlus(logZ, lw + lL)

    Xlast <- X    
  }

  return(list(logevidence=unname(logZ), posterior=posterior.samples))

}
