  dgamma_mix <- function(x, pi, alpha, beta){
    if(length(x) == 1){
      g <- length(pi)
      if(sum(pi) == 1 && min(pi) > 0 && length(alpha) == g && length(beta) == g &&
         min(alpha) > 0 && min(beta) > 0){
        aux = 0
        for(j in 1:g){aux = aux + pi[j]*dgamma(x,alpha[j], scale = beta[j])}
        return(aux)
      }
    }else{
      h <- function(x){dgamma_mix(x, pi, alpha, beta)}
      return(sapply(x, h))
    }
  }
