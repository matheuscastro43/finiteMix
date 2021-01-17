qglindley_mix <- function(p, pi, alpha, beta, gamma, lower.tail = TRUE, log.p = FALSE){
  if(length(p) == 1){
    g <- length(pi)
    if(sum(pi) == 1 && min(pi) > 0 && length(alpha) == g && length(beta) == g && 
       length(gamma) == g && min(c(alpha, beta, gamma)) > 0){
      if(log.p){
        p = exp(p)
      }
      if(!lower.tail){
        p = 1 - p
        }
      if(p < 0 || p > 1){
        warning("The probability must be in (0, 1) interval.")
        return(NaN)
      }
      else{
        if(p == 0){
          return(0)
        }
        else{
          if(p == 1){
            return(Inf)
          }
          else{
            U = aux = 100
            h = function(q){
              pglindley_mix(q, pi, alpha, beta, gamma) - p
            }
            while(aux >= 0.9*U){
              U = 2*U
              aux = uniroot(h, lower = 0, upper = U)$root
            }
            return(aux)
          }
        }
      }
    }
    else stop("The parametric space must be respected.")
    
  }else{
    j = function(p){qglindley_mix(p, pi, alpha, beta, gamma, lower.tail, log.p)}
    return(sapply(p, j))
  }
}
