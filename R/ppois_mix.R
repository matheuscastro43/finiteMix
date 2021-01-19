ppois_mix <- function(q, pi, lambda, lower.tail = TRUE, log.p = FALSE){
  if(length(q) == 1){
    g = length(pi)
    if(sum(pi) == 1 && min(c(pi, lambda)) > 0 && length(lambda) == g){
      aux = 0
      for(j in 1:g){aux = aux + pi[j]*ppois(q,lambda=lambda[j])}
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
    h = function(q){ppois_mix(q, pi, lambda, lower.tail, log.p)}
    return(sapply(q, h))
  }
}
