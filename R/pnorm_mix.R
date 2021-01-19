pnorm_mix <- function(q, pi, mean, sd, lower.tail = TRUE, log.p = FALSE){
  if(length(q) == 1){
    g = length(pi)
    if(sum(pi) == 1 && min(c(pi, sd)) > 0 && length(mean) == g && length(sd) == g){
      aux = 0
      for(j in 1:g){aux = aux + pi[j]*pnorm(q,mean=mean[j], sd = sd[j])}
      if(!lower.tail){
        aux = 1 - aux
      }
      if(!log.p){
        return(aux)
      }else{
        return(log(aux))
      }
    }else{
      stop("The parametric space must be respected.")
    }
  }else{
    h = function(q){pnorm_mix(q, pi, mean, sd, lower.tail, log.p)}
    return(sapply(q, h))
  }
}
