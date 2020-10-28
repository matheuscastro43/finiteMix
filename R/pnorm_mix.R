pnorm_mix <- function(q, pi, mean, sd, lower.tail = TRUE){
  if(length(q) == 1){
    g <- length(pi)
    if(sum(pi) == 1 && min(pi) > 0 && length(mean) == g && length(sd) == g && min(sd)>0){
      aux = 0
      for(j in 1:g){aux = aux + pi[j]*pnorm(q,mean=mean[j], sd = sd[j])}
      if(lower.tail){return(aux)
      }else{return(1 - aux)}
    }
  }else{
    h <- function(q){pnorm_mix(q, pi, mean, sd, lower.tail)}
    return(sapply(q, h))
  }
}
