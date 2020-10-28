pgamma_mix <- function(q, pi, alpha, beta, lower.tail = TRUE){
  if(length(q) == 1){
    g <- length(pi)
    if(sum(pi) == 1 && min(pi) > 0 && length(alpha) == g && length(beta) == g &&
       min(alpha) > 0 && min(beta) > 0){
      aux = 0
      for(j in 1:g){aux = aux + pi[j]*pgamma(q, alpha[j], scale = beta[j])}
      if(lower.tail){return(aux)
      }else{return(1 - aux)}
    }
  }else{
    h <- function(q){pgamma_mix(q, pi, alpha, beta, lower.tail)}
    return(sapply(q, h))
  }
}