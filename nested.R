# nessar.R - nested sampling in R
#
# Code by Matthew Hartley <Matthew.Hartley@jic.ac.uk>
# Ported from FORTRAN version by Nick Pullen and Richard Morris

library(deSolve)

initial.conditions <- c(5, 0, 15, 0, 0, 0)

pars <- c(alpha = 2,
          beta = 125)

times <- seq(0, 100, by = 0.1)

#out <- ode(initial.conditions, times, repress, pars)

debug_message <- function(message) {
  # Output message to the console

  cat(sprintf("%s\n", message))
}



evidence <- function(   prior.size,
                        prior.min,
                        prior.max,
                        n.params,
                        inital.conditions,
                        times,
                        true.data,
                        system,
                        step,
                        lower.bounds,
                        upper.bounds)
                    {
  # TODO - should prior distributions be different for different parameters?
  # TODO - this would be nicer if we just passed in the prior already
  # Create the prior
  prior.samples = replicate(prior.size, prior(prior.min, prior.max, n.params))
  llFun = llMaker(initial.conditions, times, true.data, system)
  ll.values = apply(prior.samples, 2, llFun)
  evaluated.samples = rbind(ll.values, prior.samples)

  # Order the samples by the log likelihood, worst first (at index 1)
  ordered.samples <- evaluated.samples[,order(evaluated.samples[1,])]

  # Initialise an empty matrix to hold posterior samples & their loglikelihoods
  posterior.samples <- matrix(nrow=3)
  first <- TRUE

  while ( dim(posterior.samples)[2] < 40 ) {
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
    ret <- explore(ordered.samples[2:3,selected.point], step, llMin, lower.bounds, upper.bounds)
    step <- ret$new.step
    new.point <- ret$new.values
    new.ll <- ret$new.ll
  
    # Replace the worst sample by the results of explore()
    ordered.samples[,1] <- c(new.ll, new.point)
  
    # Re-sort the samples (not very efficient, but we don't care)
    ordered.samples <- ordered.samples[,order(ordered.samples[1,])]
  }

}

MCMCStep <- function(current.values, llFun, llMin, step, lower.bounds, upper.bounds) {
  # Make single MCMC step from current.values as follows:
  # 1) Attempt step in each parameter direction
  # 2) If step is successful (ll bigger than llMin), note this
  accepted <- rep(FALSE, length(pars))
  for (i in 1:length(pars)) {
    empty <- rep(0, length(pars))
    empty[i] <- step[i]
    # Attempt step along one axis in parameter space
    candidate <- BoundedStep(current.values, empty, lower.bounds[i], upper.bounds[i])
    if(llFun(candidate) > llMin) {
      current.values <- candidate
      accepted[i] = TRUE
    }
  } 
  return(list(accepted=accepted, new.values=current.values))
}

explore <- function(current.values, step, llMin, lower.bounds, upper.bounds) {
  # Explore parameter space around supplied point
  # returns new point and new step size
  msub <- 20
  # FIXME - hardcoded the system below...
  llFun = llMaker(initial.conditions, times, true.data, repress)

  accepted <- vector()
  for(i in 1:msub) {
    ret <- MCMCStep(current.values, llFun, llMin, step, lower.bounds, upper.bounds)
    accepted <- cbind(accepted, ret$accepted)
    current.values <- ret$new.values
  }

  # DEBUG
  cat(sum(accepted[1,]), " ", sum(accepted[2,]), "\n")

  for(i in 1:length(current.values)) {
    ratio <- sum(accepted[i,]) / length(accepted[i,])
    if (ratio > 0.6) {
        step[i] <- step[i] * (1 + 2 * (ratio - 0.6) / 0.4)
    }
    if (ratio < 0.4) {
        step[i] <- step[i] / (1 + 2 * ((0.4 - ratio) / 0.4))
    }
    step[i] <- min(step[i], 0.1 * (upper.bounds[i] - lower.bounds[i]))

    #cat(sprintf("%f, %f\n", ratio, step[i]))
  }

  return(list(new.values=current.values, new.step=step, new.ll=llFun(current.values)))
}

llMaker <- function(    initial.conditions,
                        times,
                        true.data,
                        system)
{
# Return closure of logLikelihood such that it becomes a function of parameters only.
# Used so that it can be mass applied to a set of parameter choices.
    ll <- function(theta) {
        return(logLikelihood(theta, initial.conditions, times, true.data, system))
    }

    return(ll)
}

MakeStep <- function(current.values, step) {
  # Make random step in parameter space
  return(current.values + step * runif(1, -1, 1))
}

BoundedStep <- function(current.values, step, lower.bound, upper.bound) {
  # Make step in parameter space. If step exceeds bounds, instead pick point in
  # space randomly chosen using uniform distribution.

  step.attempt = MakeStep(current.values, step)
  # Generate vector of booleans that are true where we exceed bounds
  outside.bounds = (step.attempt < lower.bound) + (step.attempt > upper.bound) > 0
  # Replace values outside the bounds by randomly chosen point in parameter space
  step.attempt[which(outside.bounds)] = runif(1, lower.bound, upper.bound)

  return(step.attempt)
}

logLikelihood <- function(  theta, 
                            initial.conditions,
                            times,
                            true.data,
                            system)
{
    llh.sigma.factor <- 0.1
    range <- max(true.data) - min(true.data)
    sigma <- range * llh.sigma.factor
    guessed.data <- solveSystem(times, initial.conditions, theta, system)
    H <- -(sum((guessed.data - true.data) ** 2))/ (2 * sigma ** 2)

    return(H)
}
    

solveSystem <- function(   times, 
                            initial.conditions, 
                            pars,
                            system)
{
# Solve the model given the parameters passed in

    ret <- ode(initial.conditions, times, system, pars)

    return(ret[,4])

}

prior <- function(  xmin, xmax, n) {
# Return a vector of length n samples from the uniform distribution in the
# range (xmin, xmax)

    return( runif(n, xmin, xmax) )
}

#initialise
true.data <- read.table('true-rep-data', header=T)$gfp
