pcnorm <- function(q, pi, mean, sd, gamma, lower.tail = TRUE, log.p = FALSE){
  if(length(q) == 1){
    if(sum(pi) == 1 && min(c(pi, sd, gamma)) > 0 && length(mean) == 1 && length(sd) == 1 &&
       length(gamma) == 1 && length(pi) == 2){
      aux = pi[1]*pnorm(q, mean, sd = sd/sqrt(gamma)) + pi[2]*pnorm(q, mean, sd = sd)
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
    h = function(q){pcnorm(q, pi, mean, sd, gamma, lower.tail, log.p)}
    return(sapply(q, h))
  }
}
