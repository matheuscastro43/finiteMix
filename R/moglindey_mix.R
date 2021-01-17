moglindley_mix <- function(pi, alpha, beta, gamma){
  g <- length(pi)
  if(length(alpha) == g && length(beta) == g && length(gamma) == g && sum(pi) == 1 &&
     min(c(pi, alpha, beta, gamma)) > 0){
        f <- function(x) -dglindley_mix(x, pi, alpha, beta, gamma)
        aux = optimize(f, interval = c(0, 100))$minimum
      return(aux)
  }
  else{cat("Error.")}
}
