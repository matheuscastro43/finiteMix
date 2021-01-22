rbetar = function(n, pi, mu, phi, plot.it = TRUE, empirical = FALSE, col.pop = "red3",
                  col.empirical = "navy", ...){
  if(n == floor(n) && sum(pi) == 1 && min(c(pi, mu, phi, n)) > 0 && length(mu) == 1 && length(phi) == 1 &&
     length(pi) == 2 && mu < 1){
    
    z = rmultinom(n = n, size = 1, pi)
    aux = rowSums(z)
    modal = dbetar(optimize(function(x) dbetar(x, pi, mu, phi), interval= c(0,1), maximum = T)$maximum, pi, mu, phi)
    
    sample = c(runif(aux[1]), 
               rbeta(aux[2], mu*phi, (1-mu)*phi + 1))
    
    if(plot.it){
      d.breaks <- ceiling(nclass.Sturges(sample)*2.5)
      modal = min(c(1, max(modal, hist(sample, if(any(names(list(...)) == "breaks") == FALSE){
        breaks = d.breaks}, ...)$density)))
      hist(sample,freq = F,border = "gray48",
           main = "Sampling distribution of X",xlab = "x",
           ylab = "Density",
           xlim = c(0, 1),
           ylim = c(0, modal),
           if(any(names(list(...)) == "breaks") == FALSE){breaks = d.breaks}, ...)
      pop = function(x){dbetar(x, pi, mu, phi)}
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
      output = list(sample, pi, mu, phi, p)
      names(output) = c("sample", "pi", "mu", "phi", "plot")
    }
    else{
      output = list(sample, pi, mu, phi)
      names(output) = c("sample", "pi", "mu", "phi")
    }
    return(output)}
  else stop("The parametric space must be respected.")
}
