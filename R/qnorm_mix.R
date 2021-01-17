qnorm_mix <- function(p, pi, mean, sd, lower.tail = TRUE){
  if(length(p) == 1){
    g <- length(pi)

    if(sum(pi) == 1 && min(pi) > 0 && length(mean) == g &&
       length(sd) == g && min(sd) > 0){
      if(lower.tail == FALSE){p = 1 - p}
      if(p < 0 || p > 1){return(NaN)}
      else{
        if(p == 0){return(-Inf)}
        else{
          if(p == 1){return(Inf)}
          else{
            h <- function(q){abs(pnorm_mix(q, pi, mean, sd) - p)}
            return(optim(par = mean(mean), h, method = "L-BFGS-B")$par)
          }
        }
      }
    }

  }else{
    j <- function(p){qnorm_mix(p, pi, mean, sd, lower.tail)}
    return(sapply(p, j))
  }
}
