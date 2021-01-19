dexp_mix <- function(x, pi, rate, log = FALSE){
  if(length(x) == 1){
    g = length(pi)
    if(sum(pi) == 1 && min(c(pi, rate)) > 0 && length(rate) == g){
      aux = 0
      if(x >= 0){
        for(j in 1:g){aux = aux + pi[j]*dexp(x, rate=rate[j])}
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
    h = function(x){dexp_mix(x, pi, rate, log)}
    return(sapply(x, h))
  }
}
