dglindley <- function(x, alpha, beta, gamma){
  if(length(x) == 1){
    if(min(c(alpha, beta, gamma)) > 0){
      aux = 0
      pi = c(1/(1 + beta*gamma), beta*gamma/(1 + beta*gamma))
      return(dgamma_mix(x, pi, c(alpha, alpha + 1), rep(beta, 2)))
    }
  }else{
    h <- function(x){dglindley(x, alpha, beta, gamma)}
    return(sapply(x, h))
  }
}
