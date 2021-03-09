dcnorm <- function(x, pi, mean, sd, gamma, log = FALSE){
  if(length(x) == 1){
    pi = pi/sum(pi)
    if(length(pi) == 2 && length(mean) == 1 && length(sd) == 1 &&
       length(gamma) == 1 && min(c(pi, sd, gamma)) > 0){
      aux = pi[1]*dnorm(x, mean, sd = sd/sqrt(gamma)) + pi[2]*dnorm(x, mean, sd = sd)
      if(!log){
        return(aux)
      }else{
        return(log(aux))
      }
    }else{
      stop("The parametric space must be respected.")
    }
  }else{
    h = function(x){dcnorm(x, pi, mean, sd, gamma, log)}
    return(sapply(x, h))
  }
}
