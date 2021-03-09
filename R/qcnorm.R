qcnorm <- function(p, pi, mean, sd, gamma, lower.tail = TRUE, log.p = FALSE){
  if(length(p) == 1){
    pi = pi/sum(pi)
    if(min(c(pi, sd, gamma)) > 0 && length(mean) == 1 &&
       length(sd) == 1 && length(gamma) == 1 && length(pi) == 2){
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
              pcnorm(q, pi, mean, sd, gamma) - p
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
    j = function(p){qcnorm(p, pi, mean, sd, gamma, lower.tail, log.p)}
    return(sapply(p, j))
  }
}
