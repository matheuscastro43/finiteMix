dgamma_mix <- function(x, pi, alpha, beta, log = FALSE){
  if(length(x) == 1){
    g = length(pi)
    pi = pi/sum(pi)
    if(min(c(pi, alpha, beta)) > 0 && length(alpha) == g && 
       length(beta) == g){
      aux = 0
      if(x >= 0){
        for(j in 1:g){aux = aux + pi[j]*dgamma(x,alpha[j], scale = beta[j])}
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
    h = function(x){dgamma_mix(x, pi, alpha, beta, log)}
    return(sapply(x, h))
  }
}
