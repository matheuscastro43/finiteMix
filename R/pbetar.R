pbetar <- function(q, pi, mu, phi, lower.tail = TRUE){
  if(length(q) == 1){
    g <- length(pi)
    if(sum(pi) == 1 && min(pi) > 0 && length(pi) == 2 && length(mu) == 1 && length(phi) == 1){
      aux = pi[1] + pi[2]*pbeta(q, mu*phi, (1-mu)*phi)
      if(lower.tail){return(aux)
      }else{return(1 - aux)}
    }
  }else{
    h <- function(q){pbetar(q, pi, mu, phi, lower.tail)}
    return(sapply(q, h))
  }
}
