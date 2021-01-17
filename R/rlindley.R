rlindley = function(n, beta, plot.it = TRUE, empirical = FALSE, col.pop = "red3",
                     col.empirical = "navy", ...){
  if(n == floor(n) && n > 0 && min(beta) > 0 && is.logical(plot.it) && is.logical(empirical)
  ){
    pi = c(1/(1 + beta), beta/(1 + beta))
    z = rmultinom(n = n, size = 1, pi)
    aux = rowSums(z)
    modal <- max(dlindley(c(mogamma(1, beta), mogamma(2, beta)), beta))
    
    sample = rgamma_mix(n, pi, c(1, 2), rep(beta, 2), plot.it = FALSE)$sample
    if(plot.it){
      d.breaks <- ceiling(nclass.Sturges(sample)*2.5)
        modal = max(modal, max(density(sample)$y))
      hist(sample, freq = F, border = "gray48",
           main = "Sampling distribution of X", xlab = "x",
           ylab = "Density",
           ylim = c(0, modal),
           if(any(names(list(...)) == "breaks") == FALSE){breaks = d.breaks}, ...)
      pop = function(x){dlindley(x, beta)}
      curve(pop, col = col.pop, lwd = 3, add = T)
      if(empirical){
        lines(density(sample),col = col.empirical,lwd = 3)
        legend("topright", legend=(c("Population", "Empirical")),
               fill=c(col.pop, col.empirical),border = c(col.pop, col.empirical), bty="n")
      }
      else{
        legend("topright", legend=(c("Population")),
               fill=c(col.pop),border = c(col.pop), bty="n")
      }
      p <- recordPlot()
    }
    sample <- sort(sample)
    if(plot.it){
      output = list(sample, beta, p)
      names(output) = c("sample", "beta", "plot")
    }
    else{
      output = list(sample, beta, sample[,2])
      names(output) = c("sample", "beta")
    }
    return(output)}
  else{
    stop("Error.")
  }
}
