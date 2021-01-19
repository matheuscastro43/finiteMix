pweibull_mix <- function(q, pi, shape, scale, lower.tail = TRUE, log.p = FALSE){
  if(length(q) == 1){
    g = length(pi)
    if(sum(pi) == 1 && min(c(pi, shape, scale)) > 0 && length(shape) == g &&
       length(scale) == g){
      aux = 0
      for(j in 1:g){aux = aux + pi[j]*pweibull(q,shape=shape[j], scale = scale[j])}
      if(!lower.tail){
        aux = 1 - aux
      }
      if(!log.p){
        return(aux)
      }else{
        return(log(aux))
      }
    }else{
      stop("The parametric space must be respected.")
    }
  }
  else{
    h = function(q){pweibull_mix(q, pi, shape, scale, lower.tail, log.p)}
    return(sapply(q, h))
  }
}
