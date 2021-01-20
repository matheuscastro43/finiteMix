rglindley = function(n, alpha, beta, gamma, plot.it = TRUE, empirical = FALSE, 
                     col.pop = "red3", col.empirical = "navy", ...){
  if(n == floor(n) && min(c(alpha, beta, gamma, n)) > 0){
    pi = c(1/(1 + beta*gamma), beta*gamma/(1 + beta*gamma))
    z = rmultinom(n = n, size = 1, pi)
    aux = rowSums(z)
    modal = max(dglindley(moglindley(alpha, beta, gamma), alpha, beta, gamma))
    
    sample = rgamma_mix(n, pi, c(alpha, alpha + 1), rep(beta, 2), plot.it = FALSE)$sample
    if(plot.it){
      d.breaks <- ceiling(nclass.Sturges(sample)*2.5)
      modal = max(modal, max(density(sample)$y))
      hist(sample,freq = F,border = "gray48",
           main = "Sampling distribution of X",xlab = "x",
           ylab = "Density",
           ylim = c(0, modal),
           if(any(names(list(...)) == "breaks") == FALSE){breaks = d.breaks}, ...)
      pop = function(x){dglindley(x, alpha, beta, gamma)}
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
      output = list(sample, alpha, beta, gamma, p)
      names(output) = c("sample", "alpha", "beta", "gamma", "plot")
    }
    else{
      output = list(sample, alpha, beta, gamma)
      names(output) = c("sample", "alpha", "beta", "gamma")
    }
    return(output)}
  else stop("The parametric space must be respected.")
}
