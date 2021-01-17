qglindley <- function(p, alpha, beta, gamma, lower.tail = TRUE){
  if(length(p) == 1){

    if(min(c(alpha, beta, gamma)) > 0){
      if(lower.tail == FALSE){p = 1 - p}
      if(p < 0 || p > 1){return(NaN)}
      else{
        if(p == 0){return(0)}
        else{
          if(p == 1){return(Inf)}
          else{
            h <- function(q){abs(pglindley(q, alpha, beta, gamma) - p)}
            return(optim(par = moglindley(alpha, beta, gamma), h, method = "L-BFGS-B")$par)
          }
        }
      }
    }
    
  }else{
    j <- function(p){qglindley(p, alpha, beta, gamma, lower.tail)}
    return(sapply(p, j))
  }
}
