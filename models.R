# m = 0.1, c = 0.0003

sigmoidModel <- function(params, inputdata) {
  # Computes the output of a sigmoidal model of the inputdata
  #
  # Args:
  #  params: A vector of length four containing the parameters.
  #  inputdata: A vector of the input data to be modelled.
  #
  # Returns:
  #  The predicted data according to the model k4 + (k4 - k1) / (1 + exp(-k2 * (t - k3)))

  k1 <- params[1]
  k2 <- params[2]
  k3 <- params[3]
  k4 <- params[4]

  return(k4 + (k4 - k1) / (1 + exp(-k2 * (inputdata - k3))))
}

linearModel <- function(params, inputdata) {
  # Computes the output of a linear model of the inputdata
  #
  # Args:
  #  params: A vector of length two containing the parameters m and c.
  #  inputdata: A vector of the input data to be modelled.
  #
  # Returns:
  #  The predicted data according to the model y = m * x + c

  m <- params[1]
  c <- params[2]

  return(m * inputdata + c)
}

quadraticModel <- function(params, inputdata) {
  # Computes the output of a quadratic model of the inputdata
  #
  # Args:
  #  params: A vector of length three containing the parameters alpha, beta and gamma.
  #  inputdata: A vector of the input data to be modelled.
  #
  # Returns:
  #  The predicted data according to the model y = alpha*x^2 + beta * x + gamma

  alpha <- params[1]
  beta <- params[2]
  gamma <- params[3]

  return(alpha*inputdata*inputdata + beta*inputdata + gamma)
}

generalLogLikelihood <- function(model, params, data, llh.sigma.factor) {
  # Calculates the log likelihood using sigma factors for a given model
  #
  # Args:
  #  model: A function taking the parameters and first column of data matrix
  #         as input, and returning the predicted data.
  #  params: A vector of the parameters for the model.
  #  data: A data set containing the data to be fitted in the first column and
  #        the observed values in the second column.
  #  llh.sigma.factor: A scalar giving the sigma factor to be used
  #
  # Returns:
  #  H: The log likelihood
  
  pred <- model(params, as.matrix(data[1]))

  actual <- as.matrix(data[2])

  range <- max(actual) - min(actual)
  sigma <- range * llh.sigma.factor

  all_ll <- (actual - pred) ** 2

  H <- -sum(all_ll) / (2 * (sigma ** 2))

  return(H)

}

gradLogLikelihood <- function(model, params, data, llh.sigma.factor) {
  # IMPORTANT: THIS COULD VERY LIKELY BE WRONG!!!
  # Calculates the gradient of the log likelihood
  #
  # Args:
  #  model: A function taking the parameters and first column of data matrix
  #         as input, and returning the predicted data.
  #  params: A vector of the parameters for the model.
  #  data: A data set containing the data to be fitted in the first column and
  #        the observed values in the second column.
  #  llh.sigma.factor: A scalar giving the sigma factor to be used
  #
  # Returns:
  #  grad: The gradient

  pred <- model(params, as.matrix(data[1]))
  actual <- as.matrix(data[2])

  sigma <- 1

  grad <- (actual - pred) / sigma ** 2

  return(-sum(grad))
}
