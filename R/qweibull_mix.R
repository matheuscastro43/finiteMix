qweibull_mix <- function(p, pi, shape, scale, lower.tail = TRUE){
  if(length(p) == 1){
    g <- length(pi)
    
    if(sum(pi) == 1 && min(pi) > 0 && length(shape) == g &&
       length(scale) == g && min(c(shape, scale)) > 0){
      if(lower.tail == FALSE){p = 1 - p}
      if(p < 0 || p > 1){return(NaN)}
      else{
        if(p == 0){return(0)}
        else{
          if(p == 1){return(Inf)}
          else{
            h <- function(q){abs(pweibull_mix(q, pi, shape, scale) - p)}
            return(optim(par = mean(scale * gamma(1 + 1/shape)), h, method = "L-BFGS-B")$par)
          }
        }
      }
    }
    
  }else{
    j <- function(p){qweibull_mix(p, pi, shape, scale, lower.tail)}
    return(sapply(p, j))
  }
}
