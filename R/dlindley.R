dlindley <- function(x, beta, log = FALSE){
  if(length(x) == 1){
    if(beta > 0 && length(beta) == 1){
      aux = 0
      if(x >= 0){
        aux = (x + 1)/(beta*(beta + 1))*exp(-x/beta)
      }
      if(!log){
        return(aux)
      }else{
        return(log(aux))
      }
    }else{
      stop("The parametric space must be respected.")
    }
  }else{
    h = function(x){dlindley(x, beta, log)}
    sapply(x, h)
  }
}
