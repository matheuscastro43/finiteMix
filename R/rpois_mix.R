rpois_mix = function(n, pi, lambda, plot.it = TRUE, empirical = FALSE, col.pop = "red3",
                     col.empirical = "navy", ...){
  g = length(pi)
  pi = pi/sum(pi)
  if(n == floor(n) && min(c(pi, lambda, n)) > 0 && length(lambda) == g){
    
    z = rmultinom(n = n, size = 1, pi)
    aux = rowSums(z)
    
    sample = NULL
    for(j in 1:g){
      sample = c(sample, rpois(aux[j], lambda = lambda[j]))
    }
    if(plot.it){
      d.breaks <- ceiling(nclass.Sturges(sample)*2.5)
      modal = max(max(dpois_mix(floor(lambda), pi, lambda)),
                  hist(sample, plot = FALSE, 
                       if(any(names(list(...)) == "breaks") == FALSE){
                         breaks = d.breaks}, ...)$density)
      hist(sample,freq = F,border = "gray48",
           main = "Sampling distribution of X",xlab = "x",
           ylab = "Density", lwd = 2,
           ylim = c(0, modal),
           if(any(names(list(...)) == "breaks") == FALSE){breaks = d.breaks}, ...)
      pop = function(x){dpois_mix(x, pi, lambda)}
      x0 <- min(sample)
      x1 <- max(sample)
      par(new = T)
      plot(x0:x1, lapply(x0:x1, pop), type = "h", col = col.pop, lwd = 3,
           main = "", xlab = "",
           ylab = "", axes = F)
      if(empirical){
        lines(density(sample), col = col.empirical, lwd = 2)
        legend("topright", legend=(c("Population", "Empirical")),
               fill=c(col.pop, col.empirical),border = c(col.pop, col.empirical), bty="n")
      }
      else{
        legend("topright", legend=(c("Population")),
               fill=c(col.pop),border = c(col.pop), bty="n")
      }
      p <- recordPlot()
    }
    ord <- order(sample)
    sample <- cbind(sample, rep(1:g, aux))
    sample <- sample[ord,]
    if(plot.it){
      output = list(sample[,1], g, pi, lambda, sample[,2], p)
      names(output) = c("sample", "g", "pi", "lambda", "classification", "plot")
    }else{
      output = list(sample[,1], g, pi, lambda, sample[,2])
      names(output) = c("sample", "g", "pi", "lambda", "classification")
    }
    return(output)}
  else stop("The parametric space must be respected.")
}
