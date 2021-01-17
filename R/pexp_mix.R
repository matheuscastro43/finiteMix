pexp_mix <- function(q, pi, rate, lower.tail = TRUE){
  if(length(q) == 1){
    g <- length(pi)
    if(sum(pi) == 1 && min(pi) > 0 && length(rate) == g && min(rate) > 0){
    aux <- 0
    if(q > 0){
      for(j in 1:g){aux = aux + pi[j]* pexp(q, rate = rate[j])}
    }
    if(lower.tail){return(aux)
    }else{return(1 - aux)}
    }
  }else{
    h <- function(q){pexp_mix(q, pi, rate, lower.tail)}
    return(sapply(q, h))
  }
}
