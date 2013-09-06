baranyiModel <- function(params, inputdata) {
  # Computes the output of the Jozsef Baranyi model for bacterial growth from the inputdata
  #
  # Args:
  #  params: A vector of length four containing the parameters.
  #  inputdata: A vector of the input data to be modelled.
  #  NB: Two parameters, v and m are hardcoded to 1 for now on Jozsef's advice
  #
  # Returns:
  #  The predicted data according to the model y(t) = y0 + mu*A -(1.0/m)*log(1+((exp(m*mu*A) - 1)/(exp(m*(ymax-y0)))))
  #  where A(t) = t - lambda + log(1.0 - exp(-v*t) + exp(-v*(t-lambda)))/v

  v <- 1
  m <- 1

  y0 <- params[1]
  ymax <- params[2]
  mu <- params[3]
  lambda <- params[4]

  t <- inputdata

  A <- t - lambda + log(1.0 - exp(-v*t) + exp(-v*(t-lambda)))/v

  return(y0 + mu*A -(1.0/m)*log(1+((exp(m*mu*A) - 1)/(exp(m*(ymax-y0))))))
}

generalLogLikelihood <- function(model, params, data, weights=NULL) {
  # Calculates the log likelihood using sigma factors for a given model
  #
  # Args:
  #  model: A function taking the parameters and first column of data matrix
  #         as input, and returning the predicted data.
  #  params: A vector of the parameters for the model.
  #  data: A data set containing the data to be fitted in the first column and
  #        the observed values in the second column.
  #  weights: A vector giving the weights of each data point. If not given we should 
  #        probably integrate it out rather than set it. Less believeable points should
  # 	   take a lower weight.
  # 
  # CHECK:
  #  What I wrote for weights is true i.e. lower value.
  #
  # Returns:
  #  logLH: The log likelihood
  
  pred <- model(params, as.matrix(data[1]))

  actual <- as.matrix(data[2])

  # This needs to be done better!!! 1 may or may not be appropriate for any data set!
  if (is.null(weights)) {
    weights <- rep(length(data),1)
  }

  sigma <- weights

  all_ll <- ((actual - pred) ** 2) /  (2 * (sigma ** 2))

  logLH <- -sum(all_ll)

  return(logLH)

}
