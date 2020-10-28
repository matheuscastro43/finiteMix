dcnorm <- function(x, pi, mean, sd, gamma){
  if(length(x) == 1){
    if(sum(pi) == 1 && length(pi) == 2 && length(mean) == 1 && length(sd) == 1 &&
       length(gamma) == 1 && min(c(pi, sd, gamma)) > 0){
      return(pi[1]*dnorm(x, mean, sd = sd/sqrt(gamma)) + pi[2]*dnorm(x, mean, sd = sd))
    }
  }else{
    h <- function(x){dcnorm(x, pi, mean, sd, gamma)}
    return(sapply(x, h))
  }
}
