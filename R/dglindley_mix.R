dglindley_mix <- function(x, pi, alpha, beta, gamma){
  if(length(x) == 1){
    g <- length(pi)
    if(sum(pi) == 1 && min(pi) > 0 && length(alpha) == g && length(beta) == g &&
       length(gamma) == g && min(c(alpha, beta, gamma)) > 0){
      aux = 0
      for(j in 1:g){aux = aux + pi[j]*dglindley(x, alpha[j], beta[j], gamma[j])}
      return(aux)
    }
  }else{
    h <- function(x){dglindley_mix(x, pi, alpha, beta, gamma)}
    return(sapply(x, h))
  }
}
