pexp_mix <- function(q, pi, rate, lower.tail = TRUE, log.p = FALSE){
  if(length(q) == 1){
    g = length(pi)
    pi = pi/sum(pi)
    if(min(c(pi, rate)) > 0 && length(rate) == g){
      aux = 0
      if(q > 0){
        for(j in 1:g){aux = aux + pi[j]* pexp(q, rate = rate[j])}
      }
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
    h = function(q){pexp_mix(q, pi, rate, lower.tail, log.p)}
    return(sapply(q, h))
  }
}
