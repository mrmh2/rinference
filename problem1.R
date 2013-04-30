library(ggplot2)

bool_to_coin <- function(b) {
  # convert boolean value to H or T, TRUE = H

  if(b) return('H')
  else return('T')
}

coin_to_bool <- function(c) {
  # convert H or T to boolean value, H = TRUE

  if ( c == 'H' ) return(TRUE)
  else return(FALSE)
}

generate_sequence <- function(n, p) {
  # generate a sequence of n values, each 'H' or 'T', with probability p of
  # generating a head

  seq <- runif(n) < p

  return(sapply(seq, bool_to_coin))
}

uniform_prior <- function(n, bounds) {
  # Return n samples drawn from a uniform distribution with upper and lower
  # bounds

  lower_bound <- bounds[1]
  upper_bound <- bounds[2]

  return(runif(n, lower_bound, upper_bound))

}

coin_likelihood <- function(data, model, parameters) {
  # Evaluate the likelihood of generating the data given the model and the
  # parameters. In this case, the data is a string of 'H' or 'T' values
  # representing heads or tails, and can be summarised as the number of heads
  # obtained from n values.

  bool_seq <- sapply(data, coin_to_bool)

  N <- length(bool_seq)
  R <- sum(bool_seq)

  H <- parameters

  likelihood <- (H ** R) * (1 - H) ** (N - R)

  return(likelihood)

}

slicefun <- function(inarr, n) {
    return (inarr[1:n])
}

stuff <- function() {
  prior <- seq(0, 1, 0.01)
  data <- c('H', 'H', 'T', 'T')
  data <- generate_sequence(32, 0.25)
  posterior <- sapply(prior, coin_likelihood, data=data, model=0)
  plot(posterior, type='l')
  plot_data <- cbind.data.frame(prior, posterior)
myseq <- 2 ** seq(0, 12)
md <- sapply(myseq, slicefun, inarr=data)


}


#coin_model <- function(parameters) {
