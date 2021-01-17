pglindley_mix <- function(q, pi, alpha, beta, gamma, lower.tail = TRUE){
  if(length(q) == 1){
    g <- length(pi)
    if(sum(pi) == 1 && min(c(pi, alpha, beta, gamma)) > 0 && length(alpha) == g 
       && length(beta) == g && length(gamma) == g){
      aux = 0
      for(j in 1:g){aux = aux + pi[j]*pglindley(q, alpha[j], beta[j], gamma[j])}
      if(lower.tail){return(aux)
      }else{return(1 - aux)}
    }
  }else{
    h <- function(q){pglindley_mix(q, pi, alpha, beta, gamma, lower.tail)}
    return(sapply(q, h))
  }
}
