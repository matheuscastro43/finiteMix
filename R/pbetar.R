pbetar <- function(q, pi, mu, phi, lower.tail = TRUE, log.p = FALSE){
  if(length(q) == 1){
    g = length(pi)
    pi = pi/sum(pi)
    if(min(c(pi, mu, phi)) > 0 && length(mu) == 1 && length(phi) == 1 &&
       length(pi) == 2 && mu < 1){
      if(q == Inf){
        aux = 1
      }else{
        aux = pi[1] * q * (q>=0 && q<=1) + pi[1] * (q > 1) + pi[2]*pbeta(q, mu*phi, (1-mu)*phi + 1)
      }
      if(!lower.tail){
        aux = 1 - aux
      }
      if(!log.p){
        return(aux)
      }else{
        return(log(aux))
      }
    }else{
      stop("The parametric space must be respected.")
    }
  }else{
    h = function(q){pbetar(q, pi, mu, phi, lower.tail, log.p)}
    return(sapply(q, h))
  }
}
