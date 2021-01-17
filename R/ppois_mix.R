ppois_mix <- function(q, pi, lambda, lower.tail = TRUE){
  if(length(q) == 1){
    g = length(pi)
    if(sum(pi) == 1 && min(pi) > 0 && length(lambda) == g && min(lambda)>0){
      aux = 0
      if(q == Inf){ aux = 1
      }else{
        if(q > 0){
          for(j in 1:g){aux = aux + pi[j]*ppois(q,lambda=lambda[j])}
        }
      }
      if(lower.tail){return(aux)
      }else{return(1 - aux)}
    }
  }else{
    h <- function(q){ppois_mix(q, pi, lambda, lower.tail)}
    return(sapply(q, h))
  }
}