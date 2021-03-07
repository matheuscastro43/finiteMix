rlindley = function(n, beta, plot.it = TRUE, empirical = FALSE, col.pop = "red3",
                    col.empirical = "navy", ...){
  if(n == floor(n) && min(c(beta, n)) > 0 && length(beta) == 1){
    pi = c(1/(1 + beta), 1 - 1/(1 + beta))
    z = rmultinom(n = n, size = 1, pi)
    aux = rowSums(z)
    modal = dlindley(0, beta)
    
    sample = rgamma_mix(n, pi, c(1, 2), rep(beta, 2), plot.it = FALSE)$sample
    if(plot.it){
      d.breaks <- ceiling(nclass.Sturges(sample)*2.5)
      modal = min(c(1, max(modal, hist(sample, plot = FALSE, 
                                       if(any(names(list(...)) == "breaks") ==
                                          FALSE){
                                         breaks = d.breaks}, ...)$density)))
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
      output = list(sample, beta)
      names(output) = c("sample", "beta")
    }
    return(output)}
  else stop("The parametric space must be respected.")
}
