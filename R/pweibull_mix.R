pweibull_mix <- function(q, pi, shape, scale, lower.tail = TRUE){
  if(length(q) == 1){
    g <- length(pi)
    if(sum(pi) == 1 && min(pi) > 0 && length(shape) == g &&
       length(scale) == g && min(c(shape, scale)) > 0){
      aux = 0
      for(j in 1:g){aux = aux + pi[j]*pweibull(q,shape=shape[j], scale = scale[j])}
      if(lower.tail){return(aux)
      }else{return(1 - aux)}
    }
  }
  else{
    h <- function(q){pweibull_mix(q, pi, shape, scale, lower.tail)}
    return(sapply(q, h))
  }
}
