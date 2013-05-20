# m = 0.1, c = 0.0003

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

  llh.sigma.factor <- 10000
  range <- max(actual) - min(actual)
  sigma <- range * llh.sigma.factor

  all_ll <- (actual - pred) ** 2

  H <- -sum(all_ll) / 2 * sigma ** 2

  return(H)

}

evidence <- function() {
  
  prior.size <- 10
  prior.lower.bounds <- c(0, 0)
  prior.upper.bounds <- c(0.5, 0.001)
  #prior.samples = replicate(prior.size, prior(prior.min, prior.max, n.params))
  n <- length(prior.upper.bounds)

  for (i in 1:n) {
    }


}

data <- read.csv('data/flowering-data.csv', sep=',', header=T)

params = c(0.1, 0.0003)



