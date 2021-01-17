dlindley <- function(x, beta){
  if(beta > 0){
    if(length(x) == 1){
      return((x + 1)/(beta*(beta + 1))*exp(-x/beta))
    }
    else{
      h <- function(x){dlindley(x, beta)}
      sapply(x, h)
    }
  }
}
