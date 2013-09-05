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

makeStep <- function(current.values, step.vector) {
  # Make random step in parameter space
  #
  # Args:
  #   current.values: Vector of current parameter values.
  #   step.vector: Vector of step sizes, of the same length as the parameter
  #                values.
  #
  # Returns:
  #   New values as a vector of length number.of.parameters

  stopifnot(length(current.values) == length(step.vector))

  return(current.values + step.vector * runif(1, -1, 1))
}

makeBoundedStep <- function(current.values, step, lower.bounds, upper.bounds) {
  # BROKEN, DO NOT USE
  # Make step in parameter space. If step exceeds bounds, instead pick point in
  # space randomly chosen using uniform distribution.
  #
  # Args:
  #   current.values: Vector of parameter values
  #   
  step.attempt <- makeStep(current.values, step)
  
  # Replace value outside the bounds by randomly chosen point in parameter space
#  if((step.attempt < lower.bounds) + (step.attempt > upper.bounds) > 0) {
 #   step.attempt = runif(1, lower.bounds, upper.bounds)
  #}
  #step.attempt[which(outside.bounds)] = runif(1, lower.bounds, upper.bounds)
  return(step.attempt)
}

makeBoundedSingleStep <- function(current.value, step, lower.bound, upper.bound) {
  # Make a step of single (scalar) parameter value. If step takes us outside bounds,
  # pick point in space from uniform distribution within bounds
  #
  # Args:
  #   current.value: scalar current parameter value
  #   step: scalar step size
  #   lower.bound: scalar lower bound
  #   upper.bound: scalar upper bound
  #
  # Returns:
  #   New parameter value

  step.attempt <- current.value + step * runif(1, -1, 1)

  if (step.attempt < lower.bound | step.attempt > upper.bound) {
    step.attempt <- runif(1, lower.bound, upper.bound)
  }

  return(step.attempt)
}

makeHmcStep <- function(current.values, llFun, llMin, steps, lower.bounds, upper.bounds, gradE,llCurrent) {
  # Make random step in parameter space using HMC as given by Mackay pg 388, Algo 30.1
  # but with changes taken from http://www.cs.utoronto.ca/~radford/ham-mcmc-simple
  # Radford uses U where Mackay uses E. H = E + K. Here t(p)*p/2 = K.
  # 
  # Returns:
  #   A list containing:
  #    accepted: a placeholder currently - might want to think about adaptive step sizes in future (Neal Sec. 5.4)
  #    new.values: a vector of the new parameter values after the step
  nick <- current.values # without this my code can break!! (current.values sometimes set to NA with many problems esp.llMin not existing)
  leapfrogSteps <- 10# number of leapfrog steps in HMC
  Loop = 20#100 number of iterations
  epsilon = 0.01#The stepsize to use for the leapfrog steps # NOTE: Requires lots of manual tuning along with leapfrogsteps
  grad.current.vals <- -gradE(current.values) # gradE needs to be a fn. Look into numDeriv, or work out expression for gradient e.g. dE/dt = (x-mu)/sigma^2 ??
#  E <- -llFun(current.values) # Energy fn is -llhood
  E <- -llCurrent
  candidate <- array(,dim=c(Loop, length(current.values) +1))#array of all (if exist) accepted params and llhoods in this HMC run

  # HMC loop
  for (hoop in 1:Loop) {
    p <- rnorm(length(current.values),0,1)
    Hcurrent <- t(p)%*%p/2 + E
    new.values <- current.values
    grad.new.vals <- grad.current.vals
    for (tau in 1:leapfrogSteps) {
      p <- p - epsilon*grad.new.vals/2.0 # make half step for the momentum
      new.values <- new.values + epsilon*p # Make a full step for the position
      grad.new.vals <- -gradE( new.values ) # find new gradient
      p <- p - epsilon*grad.new.vals/2.0 # make half step for the momentum
    }
    if(any(new.values < lower.bounds) + any(new.values > upper.bounds) > 0) {
      next
    }
    Enew.llFunscore <- -llFun(new.values)
    Hnew <- t(p)%*%p/2 + Enew.llFunscore

    # Acceptance based on Metropolis criteria. Should try to simplify like Radford has it...
    # if (runif(1) < exp(-dH)) {accept = 1}
    dH <- Hnew - Hcurrent
    if (dH < 0) {
      accept = 1
    } else if (runif(1) < exp(-dH)) {
      accept = 1
    } else {
      accept = 0
    }
  
    if (accept) { 
      grad.current.vals <- grad.new.vals
      current.values <- new.values
      E <- Enew.llFunscore
      # If we accept an HMC param set lets check if it satisfies the NS likelihood criterion
      if(-E > llMin) {
        candidate[hoop,] <- c(current.values,-E)
      } 
    }
  }

  accepted <- vector()#Just put this so I don't get an error!
  # If no param sets were accepted candidate will be full of NAs, so go back to original set
  if(all(is.na(candidate))) {
    current.values <- nick
  } else {
    current.values <- candidate[which.max(candidate[,ncol(candidate)]),-ncol(candidate)]
  }
  return(list(accepted=accepted, new.values=current.values))
}

makeMcmcStep <- function(current.values, llFun, llMin, steps, lower.bounds, upper.bounds) {
  # Make single MCMC step from current.values as follows:
  # 1) Attempt step in each parameter direction
  # 2) If step is successful (ll bigger than llMin), note this
  #
  # Args:
  #   current.values: Vector of current parameter values.
  #   llFun: Function to calculate loglikelihood from llFun(param.vector).
  #   llMin: Current minimum loglikelihood.
  #   steps: Vector of step values.
  #   lower.bounds: Vector of lower bounds.
  #   upper.bounds: Vector of upper bounds.
  #
  # Returns:
  #   A list containing:
  #    accepted: a vector of boolean values describing whether each parameter's
  #              step was accepted
  #    new.values: a vector of the new parameter values after the step
  
  n <- length(current.values)
  accepted <- rep(FALSE, n)
  for (i in 1:n) {
    candidate <- current.values
    # Attempt step along one axis in parameter space
    candidate[i] <- makeBoundedSingleStep(current.values[i], steps[i],
                                          lower.bounds[i], upper.bounds[i])
    if(llFun(candidate) > llMin) {
      current.values <- candidate
      accepted[i] = TRUE
    }
  }
  return(list(accepted=accepted, new.values=current.values))
}

explore <- function(current.values, steps, llMin, llFun, lower.bounds, upper.bounds, gradFun=NULL, mcmcMethod=NULL,llCurrent=NULL) {
  # Explore parameter space around supplied point
  # returns new point and new step size
  if(is.null(mcmcMethod)){
    msub <- 20
    m <- 5

    for(k in 1:m) {
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
    }
  } else {#HMC method
    if(is.null(gradFun)){stop("ERROR: You must pass a gradient method for HMC. line 159")}
    if(is.null(llCurrent)){stop("ERROR: You must pass a selected point's loglikelihood for HMC. line 160")}
    ret <- makeHmcStep(current.values, llFun, llMin, steps, lower.bounds, upper.bounds, gradFun,llCurrent)
    current.values <- ret$new.values
    #stop("Got to line 156")
  }
    # could put refined step stuff in later - look at Radford Neal's chapter section 5.4.2.4
  return(list(new.values=current.values, new.step=steps, new.ll=llFun(current.values)))
}

prior <- function(  xmin, xmax, n) {
# Return a vector of length n samples from the uniform distribution in the
# range (xmin, xmax)

    return( runif(n, xmin, xmax) )
}

generatePriorSamples <- function(n.samples, prior.lower.bounds,
                                 prior.upper.bounds)  {
  # Generate a set of prior samples from a uniform distrubtion within specified
  # bounds.
  #
  # Args:
  #   n.samples: The number of prior samples to be generated.
  #   prior.lower.bounds: A vector of lower bounds for the prior values.
  #   prior.upper.bounds: A vector of upper bounds for the prior values.
  #
  # Returns:
  #   Matrix of prior samples, of dimension n.params x n.samples

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
                               lower.bounds, upper.bounds, gradFun, mcmcMethod) {
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
    llSelected <- ordered.samples[1,selected.point]
    # Explore around that point, updating step size as we go
    llMin <- min(ordered.samples[1,])
    ret <- explore(ordered.samples[2:nrow(ordered.samples),selected.point], 
           steps, llMin, llFun, lower.bounds, upper.bounds, gradFun, mcmcMethod, llSelected)
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

nestedSampling <- function(llFun, prior.samples, bounds, posterior.size, gradient=NULL, steps=NULL, mcmcMethod=NULL) {
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
  #   bounds:         The bounds of the parameter values, as a 2 x m matrix, 
  #                   where m is the number of parameters, the first row is 
  #                   the lower bounds and the second the upper bounds.
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

  lower.bounds <- bounds[1,]
  upper.bounds <- bounds[2,]

  # If steps is not set calculate from range of bounds
  if (is.null(steps)) {
    steps <- (upper.bounds - lower.bounds) / 100
  }
  # Generate the posterior samples
  posterior.samples <- calculatePosterior(posterior.size, ordered.samples,
    steps, llFun, lower.bounds, upper.bounds, gradient, mcmcMethod)

  # Calculate the evidence from the posterior samples
  logZ <- calculateEvidence(posterior.samples, prior.size)

  return(list(logevidence=unname(logZ), posterior=posterior.samples))

}
