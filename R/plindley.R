plindley <- function(q, beta, lower.tail = TRUE){
  if(beta > 0){
    if(length(q) == 1){
      if(q == Inf){return(1)}
      if(q <= 0){return(0)}
      aux = 1 - (beta + q + 1)/(beta + 1)*exp(-q/beta)
      if(lower.tail) return(aux)
      else return(1 - aux)
    }
    else{
      h <- function(q){plindley(q, beta, lower.tail)}
      sapply(q, h)
    }
  }
}
