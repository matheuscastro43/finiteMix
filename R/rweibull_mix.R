rweibull_mix <- function(n, pi, shape, scale, plot.it = TRUE, empirical = FALSE, col.pop = "red3",
                         col.empirical = "navy", ...){
  g <- length(pi)
  if(n == floor(n) && sum(pi) == 1 && min(c(pi, shape, scale, n)) > 0 && length(shape) == g && 
     length(scale) == g){
    
    z <- rmultinom(n, 1, pi)
    aux <- rowSums(z)
    modal <- max(dweibull_mix(scale * ((shape-1)/shape)^(1/shape), pi, shape, scale))
    
    sample <- NULL
    for(j in 1:g){
      sample <- c(sample, rweibull(aux[j], shape[j], scale[j]))
    }
    if(plot.it){
      d.breaks <- ceiling(nclass.Sturges(sample)*2.5)
      modal = min(c(1, max(modal, hist(sample, if(any(names(list(...)) == "breaks") == FALSE){
        breaks = d.breaks}, ...)$density)))
      hist(sample,freq = F,border = 1000,
           main = "Sampling distribution of X",xlab = "x",
           ylab = "Density",
           ylim = c(0, modal),
           if(any(names(list(...)) == "breaks") == FALSE){breaks = d.breaks}, ...)
      pop = function(x){dweibull_mix(x, pi, shape, scale)}
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
      output = list(sample[,1], g, pi, shape, scale, sample[,2], p)
      names(output) = c("sample", "g", "pi", "a", "b", "classification", "plot")
    }
    else{
      output = list(sample[,1], g, pi, shape, scale, sample[,2])
      names(output) = c("sample", "g", "pi", "a", "b", "classification")
    }
    return(output)}
  else stop("The parametric space must be respected.")
}
