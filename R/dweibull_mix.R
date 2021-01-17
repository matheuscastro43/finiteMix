dweibull_mix <- function(x, pi, shape, scale){
  if(length(x) == 1){
    g <- length(pi)
    if(sum(pi) == 1 && min(pi) > 0 && length(shape) == g && length(scale) == g &&
       min(c(shape, scale)) > 0){
      aux = 0
      for(j in 1:g){aux = aux + pi[j]*dweibull(x, shape[j], scale[j])}
      return(aux)
    }
  }
  else{
    h <- function(x){dweibull_mix(x, pi, shape, scale)}
    return(sapply(x, h))
  }
}

