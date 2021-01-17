dnorm_mix <- function(x, pi, mean, sd){
  if(length(x) == 1){
    g <- length(pi)
    if(sum(pi) == 1 && min(pi) > 0 && length(mean) == g && length(sd) == g && min(sd)>0){
      aux = 0
      for(j in 1:g){aux = aux + pi[j]*dnorm(x,mean=mean[j], sd = sd[j])}
      return(aux)
    }
  }else{
    h <- function(x){dnorm_mix(x, pi, mean, sd)}
    return(sapply(x, h))
  }
}
