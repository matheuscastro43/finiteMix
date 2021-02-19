pglindley <- function(q, alpha, beta, gamma, lower.tail = TRUE, log.p = FALSE){
  if(length(q) == 1){
    if(length(alpha) == 1 && length(beta) == 1 && length(gamma) == 1){
      if(min(c(alpha, beta, gamma)) > 0 && length(alpha) == 1 && 
         length(beta) == 1 && length(gamma) == 1){
        pi = c(1/(1 + beta*gamma), beta*gamma/(1 + beta*gamma))
        return(pgamma_mix(q, pi, c(alpha, alpha + 1), rep(beta, 2), lower.tail, log.p))
      }else{
        stop("The parametric space must be respected.")
      }
    }else{
      k <- function(q, alpha, beta, gamma) pglindley(q, alpha, beta, gamma,
                                                     lower.tail, log.p)
      nt = max(length(alpha), length(beta), length(gamma))
      mapply(k, q, rep(alpha, length = nt), rep(beta, length = nt), 
             rep(gamma, length = nt))
    }
  }else{
    h = function(q){pglindley(q, alpha, beta, gamma, lower.tail, log.p)}
    return(sapply(q, h))
  }
}
