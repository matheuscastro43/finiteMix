moglindley <- function(alpha, beta, gamma){
  if((length(alpha) == length(beta)) && (length(alpha) == length(gamma)) &&
     min(c(alpha, beta, gamma)) > 0){
    g <- length(alpha)
    if(g == 1){
        f <- function(x) -dglindley(x, alpha, beta, gamma)
        aux = (optimize(f, interval = c(0, 100))$minimum)
      return(aux)
    }
    else{
      h <- function(a, b, c){moglindley(a, b, c)}
      mapply(h, alpha, beta, gamma)
    }
  }
  else{cat("Error.")}
}

