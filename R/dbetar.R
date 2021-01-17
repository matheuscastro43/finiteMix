dbetar <- function(x, pi, mu, phi){
  if(length(x) == 1){
    if(sum(pi) == 1 && min(pi) > 0 && length(mu) == 1 && length(phi) == 1 &&
       length(pi) == 2){
      aux = 0
      return(pi[1] + pi[2]*dbeta(x, mu*phi, (1-mu)*phi + 1))
    }
  }else{
    h <- function(x){dbetar(x, pi, mu, phi)}
    return(sapply(x, h))
  }
}