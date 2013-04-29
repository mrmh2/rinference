bool_to_coin <- function(b) {
  # convert boolean value to H or T, TRUE = H

  if(b) return('H')
  else return('T')
}

generate_sequence <- function(n, p) {
  # generate a sequence of n values, each 'H' or 'T', with probability p of
  # generating a head

  seq <- runif(n) < p

  return(sapply(seq, bool_to_coin))
}





