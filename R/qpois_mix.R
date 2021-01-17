qpois_mix <- function(p, pi, lambda, lower.tail = TRUE){
  if(length(p) == 1){
    g = length(pi)
    
    if(sum(pi) == 1 && min(pi) > 0 && length(lambda) == g && min(lambda)>0){
      if(lower.tail == FALSE){p = (1-p)}
      
      if(p < 0 || p > 1){return(NaN)
      }else{
        if(p == 1){return(Inf)
        }else{
          if(p == 0){return(0)
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
  }else{
    h = function(x){qpois_mix(x, pi, lambda, lower.tail)}
    return(sapply(p, h))
  }
}
