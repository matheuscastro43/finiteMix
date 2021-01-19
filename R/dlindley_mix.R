dlindley_mix <- function(x, pi, beta, log = FALSE){
  if(length(x) == 1){
    g = length(pi)
    if(sum(pi) == 1 && min(c(pi, beta)) > 0 && length(beta) == g){
      aux = 0
      if(x >= 0){
        for(j in 1:g){aux = aux + pi[j]*dlindley(x, beta = beta[j])}
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
    h = function(x){dlindley_mix(x, pi, beta, log)}
    return(sapply(x, h))
  }
}
