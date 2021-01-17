pglindley <- function(q, alpha, beta, gamma, lower.tail = TRUE){
  if(length(q) == 1){
    if(min(c(alpha, beta, gamma)) > 0){
      pi = c(1/(1 + beta*gamma), beta*gamma/(1 + beta*gamma))
      aux = pgamma_mix(q, pi, c(alpha, alpha + 1), rep(beta, 2))
      if(lower.tail){return(aux)
      }else{return(1 - aux)}
    }
  }else{
    h <- function(q){pglindley(q, alpha, beta, gamma, lower.tail)}
    return(sapply(q, h))
  }
}
