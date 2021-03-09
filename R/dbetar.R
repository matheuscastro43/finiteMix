dbetar <- function(x, pi, mu, phi, log = FALSE){
  if(length(x) == 1){
    pi = pi/sum(pi)
    if(min(c(pi, mu, phi)) > 0 && length(mu) == 1 && length(phi) == 1 &&
       length(pi) == 2 && mu < 1){
      aux = 0
      if(x >= 0 && x <= 1){
        aux = pi[1] + pi[2]*dbeta(x, mu*phi, (1-mu)*phi + 1)
      }
      if(!log){
        return(aux)
      }
      else{
        return(log(aux))
      }
    }
    stop("The parametric space must be respected.")
  }else{
    h = function(x){dbetar(x, pi, mu, phi, log)}
    return(sapply(x, h))
  }
}
