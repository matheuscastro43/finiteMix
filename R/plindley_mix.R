plindley_mix <- function(q, pi, beta, lower.tail = TRUE, log.p = FALSE){
  if(length(q) == 1){
    g = length(pi)
    pi = pi/sum(pi)
    if(min(c(pi, beta)) > 0 && length(beta) == g){
      aux = 0
      for(j in 1:g){aux = aux + pi[j]*plindley(q, beta = beta[j])}
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
  }else{
    h = function(q){plindley_mix(q, pi, beta, lower.tail, log.p)}
    return(sapply(q, h))
  }
}
