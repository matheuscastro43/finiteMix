qlindley <- function(p, beta, lower.tail = TRUE){
  if(length(p) == 1){
    if(beta > 0){
      if(!lower.tail){p = 1 - p}
      if(p < 0 || p > 1){return(NaN)}
      else{
        if(p == 0){return(0)}
        else{
          if(p == 1){return(Inf)}
          else{
            h <- function(q){abs(plindley(q, beta) - p)}
            return(optim(par = molindley(beta), h, method = "L-BFGS-B")$par)
          }
        }
      }
    }
    
  }else{
    j <- function(p){qlindley(p, beta, lower.tail)}
    return(sapply(p, j))
  }
}
