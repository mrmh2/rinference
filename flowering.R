# m = 0.1, c = 0.0003

# k1 / (1 + exp(-k2 * (t - k3)))

logplus <- function(a, b) {
    if (a > b) {
        sum = a + log(1+exp(b-a))
    } else {
        sum = b + log(1+exp(a-b))
    }
    return(sum)
}

MakeStep <- function(current.values, step) {
  # Make random step in parameter space
  return(current.values + step * runif(1, -1, 1))
}

BoundedStep <- function(current.values, step, lower.bounds, upper.bounds) {
  # Make step in parameter space. If step exceeds bounds, instead pick point in
  # space randomly chosen using uniform distribution.

  step.attempt <- MakeStep(current.values, step)
  # Generate vector of booleans that are true where we exceed bounds
  outside.bounds <- (step.attempt < lower.bounds) + (step.attempt > upper.bounds) > 0
  # Replace values outside the bounds by randomly chosen point in parameter space
  step.attempt[which(outside.bounds)] = runif(1, lower.bounds, upper.bounds)

  return(step.attempt)
}

MCMCStep <- function(current.values, llFun, llMin, steps, lower.bounds, upper.bounds) {
  # Make single MCMC step from current.values as follows:
  # 1) Attempt step in each parameter direction
  # 2) If step is successful (ll bigger than llMin), note this

  n <- length(current.values)
  accepted <- rep(FALSE, n)
  for (i in 1:n) {
    empty <- rep(0, n)
    empty[i] <- steps[i]
    # Attempt step along one axis in parameter space
    candidate <- BoundedStep(current.values, empty, lower.bounds, upper.bounds)
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
    ret <- MCMCStep(current.values, llFun, llMin, steps, lower.bounds, upper.bounds)
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

ftLogLikelihood <- function(model, params, data) {
  m <- params[1]
  c <- params[2]
  sd <- 1

  pred <- m * as.matrix(data[1]) + c
  actual <- as.matrix(data[2])

  #slike <- dnorm(as.matrix(data[2]), mean=pred, sd=sd, log=T)

  llh.sigma.factor <- 0.1
  range <- max(actual) - min(actual)
  sigma <- range * llh.sigma.factor

  all_ll <- (actual - pred) ** 2

  H <- -sum(all_ll) / (2 * (sigma ** 2))

  return(H)

}

GenerateAPrior <- function(prior.lower.bounds, prior.upper.bounds) {

  n <- length(prior.upper.bounds)
  prior <- numeric(n)

  for (i in 1:n) {
    prior[i] <- runif(1, prior.lower.bounds[i], prior.upper.bounds[i])
  }

  return(prior)
}

CalculateEvidence <- function(llFun, posterior.size, prior.samples) {
  
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

  logZ <- -1e300
  Xlast <- 1
  for (i in 1:posterior.size) {
    lL <- posterior.samples[1,i]
    X <- exp(-i / prior.size)
    lw <- log(Xlast - X)
    #cat(lw, ' ', lL, ' ', logZ + lw + lL, '\n')

    logZ <- logplus(logZ, lw + lL)

    Xlast <- X    
  }

  return(list(logevidence=logZ, posterior=posterior.samples))

}

data <- read.csv('data/flowering-data.csv', sep=',', header=T)

params = c(0.1, 0.0003)

FTllFun <- function(params) { 
  return(ftLogLikelihood(NULL, params, data)) 
}

prior.size <- 100
prior.lower.bounds <- c(0, 0)
prior.upper.bounds <- c(0.5, 0.001)
prior.samples <- replicate(prior.size, GenerateAPrior(prior.lower.bounds, prior.upper.bounds))
 
