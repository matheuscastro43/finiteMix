mogamma <- function(alpha, beta){
  if(length(alpha) == length(beta)){
    g <- length(alpha)
    if(g == 1){
      aux = 0
      for(j in c(1:g)){
        if(alpha[j] >= 1) {aux[j] = ((alpha[j] - 1)/beta)}
        else{
          f <- function(x) -dgamma(x, alpha[j], scale = beta[j])
          aux[j] = (optimize(f, interval = c(0, 100))$minimum)
        }
      }
      return(aux)
    }
    else{
      h <- function(a, b){mogamma(a, b)}
      mapply(h, alpha, beta)
    }
  }
  else{cat("Error.")}
}
