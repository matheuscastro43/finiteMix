qexp_mix <- function(p, pi, rate, lower.tail = TRUE){
  if(length(p) == 1){
    g <- length(pi)

    if(sum(pi) == 1 && min(pi) > 0 && length(rate) == g && min(rate) > 0){
      if(lower.tail == FALSE){p = (1-p)}

      if(p < 0 || p > 1){return(NaN)}
      else{
        if(p == 0){return(0)}
        else{
          if(p == 1){return(Inf)}
          else{
            h <- function(q){abs(pexp_mix(q, pi, rate) - p)}
            aux <- (optim(par = 1/min(rate), h, method = "L-BFGS-B"))$par
            if(aux > 0){return(aux)
            }else{return(0)}
          }
        }
      }
    }
  }else{
    j <- function(p){qexp_mix(p, pi, rate, lower.tail)}
    return(sapply(p, j))
  }
}
