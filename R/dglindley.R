dglindley <- function(x, alpha, beta, gamma, log = FALSE){
  if(length(x) == 1){
    if(length(alpha) == 1 && length(beta) == 1 && length(gamma) == 1){
      if(min(c(alpha, beta, gamma)) > 0 && length(alpha) == 1 && 
         length(beta) == 1 && length(gamma) == 1){
        aux = 0
        if(x >= 0){
          pi = c(1/(1 + beta*gamma), 1 - 1/(1 + beta*gamma))
          aux = dgamma_mix(x, pi, c(alpha, alpha + 1), rep(beta, 2))
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
      k <- function(x, alpha, beta, gamma) dglindley(x, alpha, beta, gamma,
                                                     log = log)
      nt = max(length(alpha), length(beta), length(gamma))
      mapply(k, x, rep(alpha, length = nt), rep(beta, length = nt), 
             rep(gamma, length = nt))
    }
  }else{
    h = function(x){dglindley(x, alpha, beta, gamma, log)}
    return(sapply(x, h))
  }
}
