qglindley_mix <- function(p, pi, alpha, beta, gamma, lower.tail = TRUE){
  if(length(p) == 1){
    g <- length(pi)
    
    if(sum(pi) == 1 && min(pi) > 0 && length(alpha) == g && length(beta) == g && 
       length(gamma) == g && min(c(alpha, beta, gamma)) > 0){
      if(lower.tail == FALSE){p = 1 - p}
      if(p < 0 || p > 1){return(NaN)}
      else{
        if(p == 0){return(0)}
        else{
          if(p == 1){return(Inf)}
          else{
            h <- function(q){abs(pglindley_mix(q, pi, alpha, beta, gamma) - p)}
            return(optim(par = moglindley_mix(pi, alpha, beta, gamma), h, method = "L-BFGS-B")$par)
          }
        }
      }
    }
    
  }else{
    j <- function(p){qglindley_mix(p, pi, alpha, beta, gamma, lower.tail)}
    return(mapply(j, p))
  }
}
