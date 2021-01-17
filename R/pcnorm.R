pcnorm <- function(q, pi, mean, sd, gamma, lower.tail = TRUE){
  if(length(q) == 1){
    if(sum(pi) == 1 && min(c(pi, sd, gamma)) > 0 && length(mean) == 1 && length(sd) == 1 &&
       length(gamma) == 1 && length(pi) == 2){
      aux = pi[1]*pnorm(q, mean, sd = sd/sqrt(gamma)) + pi[2]*pnorm(q, mean, sd = sd)
      if(lower.tail){return(aux)
      }else{return(1 - aux)}
    }
  }else{
    h <- function(q){pcnorm(q, pi, mean, sd, gamma, lower.tail)}
    return(sapply(q, h))
  }
}
