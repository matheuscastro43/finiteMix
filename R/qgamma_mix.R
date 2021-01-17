qgamma_mix <- function(p, pi, alpha, beta, lower.tail = TRUE){
  if(length(p) == 1){
    g <- length(pi)
    
    if(sum(pi) == 1 && min(pi) > 0 && length(alpha) == g &&
       length(beta) == g && min(alpha) > 0 && min(beta) > 0){
      if(lower.tail == FALSE){p = 1 - p}
      if(p < 0 || p > 1){return(NaN)}
      else{
        if(p == 0){return(0)}
        else{
          if(p == 1){return(Inf)}
          else{
            h <- function(q){abs(pgamma_mix(q, pi, alpha, beta) - p)}
            return(optim(par = 0, h, method = "L-BFGS-B")$par)
          }
        }
      }
    }
    
  }else{
    j <- function(p){qgamma_mix(p, pi, alpha, beta, lower.tail)}
    return(sapply(p, j))
  }
}