qnorm_mix <- function(p, pi, mean, sd, lower.tail = TRUE, log.p = FALSE){
  if(length(p) == 1){
    g <- length(pi)
    pi = pi/sum(pi)
    if(min(c(pi, sd)) > 0 && length(mean) == g && 
       length(sd) == g){
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
          return(-Inf)
        }
        else{
          if(p == 1){
            return(Inf)
          }
          else{
            U = aux = 100
            h = function(q){
              pnorm_mix(q, pi, mean, sd) - p
            }
            while(abs(aux) >= 0.9*U){
              U = 2*U
              aux = uniroot(h, lower = -U, upper = U)$root
            }
            return(aux)
          }
        }
      }
    }
    else stop("The parametric space must be respected.")
    
  }else{
    j <- function(p){qnorm_mix(p, pi, mean, sd, lower.tail)}
    return(sapply(p, j))
  }
}
