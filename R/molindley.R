molindley <- function(beta){
  if(beta > 1)
    return(beta - 1)
  if(beta > 0)
    return(0)
  cat("Error.")
}
