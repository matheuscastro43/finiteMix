dpois_mix <- function(x, pi, lambda, log = FALSE){
  if(length(x) == 1){
    g = length(pi)
    pi = pi/sum(pi)
    if(min(c(pi, lambda)) > 0 && length(lambda) == g){
      aux = 0
      if(x == floor(x)){
        for(j in 1:g){aux = aux + pi[j]*dpois(x,lambda=lambda[j])}
      }
      if(!log){
        return(aux)
      }else{
        return(log(aux))
      }
    }else{
      stop("The parametric space must be respected.")
    }
  }else{
    h = function(x){dpois_mix(x, pi, lambda, log)}
    return(sapply(x, h))
  }
}
