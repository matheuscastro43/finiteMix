dweibull_mix <- function(x, pi, shape, scale, log = FALSE){
  if(length(x) == 1){
    g = length(pi)
    if(sum(pi) == 1 && min(c(pi, shape, scale)) > 0 && length(shape) == g && 
       length(scale) == g){
      aux = 0
      for(j in 1:g){aux = aux + pi[j]*dweibull(x, shape[j], scale[j])}
      if(!log){
        return(aux)
      }else{
        return(log(aux))
      }
    }else{
      stop("The parametric space must be respected.")
    }
  }
  else{
    h = function(x){dweibull_mix(x, pi, shape, scale, log)}
    return(sapply(x, h))
  }
}
