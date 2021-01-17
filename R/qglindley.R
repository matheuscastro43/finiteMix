qglindley <- function (p, alpha, beta, gamma, lower.tail = TRUE){
  if (length(p) == 1) {
    if (min(c(alpha, beta, gamma)) > 0) {
      if (lower.tail == FALSE) {
        p = 1 - p
      }
      if (p < 0 || p > 1) {
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
  }
  else {
    j = function(p) {
      qglindley(p, alpha, beta, gamma, lower.tail)
    }
    return(sapply(p, j))
  }
}
