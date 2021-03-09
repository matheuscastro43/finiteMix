dglindley_mix <- function(x, pi, alpha, beta, gamma, log = FALSE){
  if(length(x) == 1){
    g = length(pi)
    pi = pi/sum(pi)
    if(min(c(pi, alpha, beta, gamma)) > 0 
       && length(alpha) == g && length(beta) == g && length(gamma) == g){
      aux = 0
      if(x >= 0){
        for(j in 1:g){aux = aux + pi[j]*dglindley(x, alpha[j], beta[j], gamma[j])}
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
    h = function(x){dglindley_mix(x, pi, alpha, beta, gamma, log)}
    return(sapply(x, h))
  }
}
