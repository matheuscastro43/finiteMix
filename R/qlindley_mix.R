qlindley_mix <- function(p, pi, beta, lower.tail = TRUE){
  if(length(p) == 1){
    g <- length(pi)
    
    if(sum(pi) == 1 && min(pi) > 0 && length(beta) == g && min(beta) > 0){
      if(!lower.tail){p = 1 - p}
      if(p < 0 || p > 1){return(NaN)}
      else{
        if(p == 0){return(0)}
        else{
          if(p == 1){return(Inf)}
          else{
            h <- function(q){abs(plindley_mix(q, pi, beta) - p)}
            return(optim(par = molindley_mix(pi, beta), h, method = "L-BFGS-B")$par)
          }
        }
      }
    }
    
  }else{
    j <- function(p){qlindley_mix(p, pi, beta, lower.tail)}
    return(sapply(p, j))
  }
}
