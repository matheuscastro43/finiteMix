molindley_mix <- function(pi, beta){
  if(beta > 0 && (length(pi) == length(beta))){
    modes <- sapply(beta, molindley)
    return(modes[dlindley_mix(modes, pi, beta) == max(dlindley_mix(modes, pi, beta))])
  }
  cat("Error.")
}
