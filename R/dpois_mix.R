dpois_mix <- function(x, pi, lambda){
  if(length(x) == 1){
    g = length(pi)
    if(sum(pi) == 1 && min(pi) > 0 && length(lambda) == g && min(lambda)>0){
      aux = 0
      if(x == floor(x)){
        for(j in 1:g){aux = aux + pi[j]*dpois(x,lambda=lambda[j])}
        return(aux)
      }
      else(return(0))
    }
  }else{
    h = function(x){dpois_mix(x, pi, lambda)}
    return(sapply(x, h))
  }
}
