dlindley_mix <- function(x, pi, beta){
  if(length(x) == 1){
    g <- length(pi)
    if(sum(pi) == 1 && min(pi) > 0 && length(beta) == g && min(beta) > 0){
      aux = 0
      for(j in 1:g){aux = aux + pi[j]*dlindley(x, beta = beta[j])}
      return(aux)
    }
  }else{
    h <- function(x){dlindley_mix(x, pi, beta)}
    return(sapply(x, h))
  }
}
