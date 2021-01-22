qglindley <- function (p, alpha, beta, gamma, lower.tail = TRUE, log.p = FALSE){
  if (length(p) == 1) {
    if (min(c(alpha, beta, gamma)) > 0 && length(alpha) == 1 && 
        length(beta) == 1 && length(gamma) == 1){
      if(log.p){
        p <- exp(p)
        }
      if (!lower.tail) {
        p = 1 - p
      }
      if (p < 0 || p > 1) {
        warning("The probability must be in (0, 1) interval.")
        return(NaN)
      }
      else {
        if (p == 0) {
          return(0)
        }
        else {
          if (p == 1) {
            return(Inf)
          }
          else {
            U = aux = 100
            h = function(q){
              pglindley(q, alpha, beta, gamma) - p
            }
            while(aux >= 0.9*U){
              U = 2*U
              aux = uniroot(h, lower = 0, upper = U)$root
            }
            return(aux)
          }
        }
      }
    }
    else stop("All parameters must be positive.")
  }
  else {
    j = function(p) {
      qglindley(p, alpha, beta, gamma, lower.tail, log.p)
    }
    return(sapply(p, j))
  }
}
