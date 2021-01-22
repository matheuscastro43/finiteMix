plindley <- function(q, beta, lower.tail = TRUE, log.p = FALSE){
  if(length(q) == 1){
    if(beta > 0 && length(beta) == 1){
      if(q != Inf){
        aux = (1 - (beta + q + 1)/(beta + 1)*exp(-q/beta))*(q >= 0)
      }else{
        aux = 1
      }
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
    h = function(q){plindley(q, beta, lower.tail, log.p)}
    sapply(q, h)
  }
}
