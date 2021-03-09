rcnorm = function(n, pi, mean, sd, gamma, plot.it = TRUE, empirical = FALSE, col.pop = "red3",
                  col.empirical = "navy", ...){
  pi = pi/sum(pi)
  if(n == floor(n) && length(pi) == 2 && length(mean) == 1 && length(sd) == 1 &&
     length(gamma) == 1 && min(c(pi, sd, gamma, n)) > 0){
    
    z = rmultinom(n = n, size = 1, pi)
    aux = rowSums(z)
    modal = max(dcnorm(mean, pi, mean, sd, gamma))
    
    sample = c(rnorm(aux[1], mean, sd = sd/sqrt(gamma)),
               rnorm(aux[2], mean, sd = sd))
    
    if(plot.it){
      d.breaks <- ceiling(nclass.Sturges(sample)*2.5)
      modal = min(c(1, max(modal, hist(sample, plot = FALSE,
                                       if(any(names(list(...)) == "breaks") == 
                                          FALSE){
                                         breaks = d.breaks}, ...)$density)))
      hist(sample,freq = F,border = "gray48",
           main = "Sampling distribution of X",xlab = "x",
           ylab = "Density",
           ylim = c(0, modal),
           if(any(names(list(...)) == "breaks") == FALSE){breaks = d.breaks}, ...)
      pop = function(x){dcnorm(x, pi, mean, sd, gamma)}
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
      output = list(sample, pi, mean, sd, gamma, p)
      names(output) = c("sample", "pi", "mu", "sigma", "gamma", "plot")
    }
    else{
      output = list(sample, pi, mean, sd, gamma)
      names(output) = c("sample", "pi", "mu", "sigma", "gamma")
    }
    return(output)}
  else stop("The parametric space must be respected.")
}
