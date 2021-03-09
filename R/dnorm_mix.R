dnorm_mix <- function(x, pi, mean, sd, log = FALSE){
  if(length(x) == 1){
    g = length(pi)
    pi = pi/sum(pi)
    if( min(c(pi, sd)) > 0 && length(mean) == g && length(sd) == g){
      aux = 0
      for(j in 1:g){aux = aux + pi[j]*dnorm(x,mean=mean[j], sd = sd[j])}
      if(!log){
        return(aux)
      }else{
        return(log(aux))
      }
    }else{
      stop("The parametric space must be respected.")
    }
  }else{
    h = function(x){dnorm_mix(x, pi, mean, sd, log)}
    return(sapply(x, h))
  }
}
