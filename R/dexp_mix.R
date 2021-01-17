dexp_mix <- function(x, pi, rate){
  if(length(x) == 1){
    g <- length(pi)
    if(g == floor(g) && g>0 && length(pi) == g && sum(pi) == 1 && min(pi) > 0
       && length(rate) == g && min(rate)>0){
      aux = 0
      if(x <= 0){return(aux)}
      for(j in 1:g){aux = aux + pi[j]*dexp(x, rate=rate[j])}
      return(aux)
    }
  }else{
    h <- function(x){dexp_mix(x, pi, rate)}
    return(sapply(x, h))
  }
}
