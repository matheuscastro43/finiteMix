qpois_mix <- function(p, pi, lambda, lower.tail = TRUE, log.p = FALSE){
  if(length(p) == 1){
    g = length(pi)
    pi = pi/sum(pi)
    if(min(c(pi, lambda)) > 0 && length(lambda) == g){
      if(log.p){
        p = exp(p)
      }
      if(!lower.tail){
        p = (1-p)
      }
      if(p < 0 || p > 1){
        warning("The probability must be in (0, 1) interval.")
        return(NaN)
      }else{
        if(p == 1){
          return(Inf)
        }else{
          if(p == 0){
            return(0)
          }else{
            aux = 0
            while(T){
              if(ppois_mix(aux, pi, lambda) > p){
                break
              }else{
                aux = aux + 1
              }
            }
            return(aux)
          }
        }
      }
    }
    else stop("The parametric space must be respected.")
  }else{
    h = function(x){qpois_mix(x, pi, lambda, lower.tail, log.p)}
    return(sapply(p, h))
  }
}
