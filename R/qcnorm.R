qcnorm <- function(p, pi, mean, sd, gamma, lower.tail = TRUE){
  if(length(p) == 1){
    
    if(sum(pi) == 1 && min(c(pi, sd, gamma)) > 0 && length(mean) == 1 &&
       length(sd) == 1 && length(gamma) == 1 && length(pi) == 2){
      if(lower.tail == FALSE){p = 1 - p}
      if(p < 0 || p > 1){return(NaN)}
      else{
        if(p == 0){return(-Inf)}
        else{
          if(p == 1){return(Inf)}
          else{
            h <- function(q){abs(pcnorm(q, pi, mean, sd, gamma) - p)}
            return(optim(par = mean, h, method = "L-BFGS-B")$par)
          }
        }
      }
    }
    
  }else{
    j <- function(p){qcnorm(p, pi, mean, sd, gamma, lower.tail)}
    return(sapply(p, j))
  }
}
