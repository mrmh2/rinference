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

}

coin_likelihood <- function(data, model, parameters) {
  # Evaluate the likelihood of generating the data given the model and the
  # parameters. In this case, the data is a string of 'H' or 'T' values
  # representing heads or tails, and can be summarised as the number of heads
  # obtained from n values.

  bool_seq <- sapply(data, coin_to_bool)

  N <- length(bool_seq)
  R <- sum(bool_seq)

  H <- 0.5

  likelihood <- (H ** R) * (1 - H) ** (N - R)

  return(likelihood)

}

#coin_model <- function(parameters) {
