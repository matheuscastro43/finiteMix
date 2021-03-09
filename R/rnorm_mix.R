rnorm_mix = function(n, pi, mean, sd, plot.it = TRUE, empirical = FALSE, col.pop = "red3",
                     col.empirical = "navy", ...){
  g <- length(pi)
  pi = pi/sum(pi)
  if(n == floor(n) && min(c(pi, sd, n)) > 0 && length(mean) == g && length(sd) == g){
    
    z = rmultinom(n = n, size = 1, pi)
    aux = rowSums(z)
    modal <- max(dnorm_mix(mean, pi, mean, sd))
    
    sample = NULL
    for(j in 1:g){
      sample = c(sample, rnorm(aux[j], mean = mean[j], sd = sd[j]))
    }
    if(plot.it){
      d.breaks <- ceiling(nclass.Sturges(sample)*2.5)
      modal = min(c(1, max(modal, hist(sample, plot = FALSE, 
                                       if(any(names(list(...)) == 
                                              "breaks") == FALSE){
                                         breaks = d.breaks}, ...)$density)))
      hist(sample,freq = F,border = "gray48",
           main = "Sampling distribution of X",xlab = "x",
           ylab = "Density",
           ylim = c(0, modal),
           if(any(names(list(...)) == "breaks") == FALSE){breaks = d.breaks}, ...)
      pop = function(x){dnorm_mix(x, pi, mean, sd)}
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
    ord <- order(sample)
    sample <- cbind(sample, rep(1:g, aux))
    sample <- sample[ord,]
    if(plot.it){
      output = list(sample[,1], g, pi, mean, sd, sample[,2], p)
      names(output) = c("sample", "g", "pi", "mu", "sigma", "classification", "plot")
    }
    else{
      output = list(sample[,1], g, pi, mean, sd, sample[,2])
      names(output) = c("sample", "g", "pi", "mu", "sigma", "classification")
    }
    return(output)
  }else stop("The parametric space must be respected.")
}
