plindley_mix <- function(q, pi, beta, lower.tail = TRUE){
  if(length(q) == 1){
    g <- length(pi)
    if(sum(pi) == 1 && min(pi) > 0 && length(beta) == g && min(beta) > 0){
      if(q == Inf){return(1)}
      if(q <= 0){return(0)}
      aux = 0
      for(j in 1:g){aux = aux + pi[j]*plindley(q, beta = beta[j])}
      if(lower.tail){return(aux)
      }else{return(1 - aux)}
    }
  }else{
    h <- function(q){plindley_mix(q, pi, beta, lower.tail)}
    return(sapply(q, h))
  }
}
