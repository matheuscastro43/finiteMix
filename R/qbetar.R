qbetar <- function(p, pi, mu, phi, lower.tail = TRUE, log.p = FALSE){
  if(length(p) == 1){
    pi = pi/sum(pi)
    if(min(c(pi, mu, phi)) > 0 && length(mu) == 1 && length(phi) == 1 &&
       length(pi) == 2 && mu < 1){
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
        if(p == 0 || p == 1){
          return(p)
        }
        else{
          h = function(q){
            pbetar(q, pi, mu, phi) - p
          }
          return(uniroot(h, lower = 0, upper = 1)$root)
        }
      }
    }else stop("The parametric space must be respected.")
  }else{
    j = function(p){qbetar(p, pi, mu, phi, lower.tail, log.p)}
    return(sapply(p, j))
  }
}
