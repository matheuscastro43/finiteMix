rglindley_mix = function(n, pi, alpha, beta, gamma, plot.it = TRUE, empirical = FALSE, 
                         col.pop = "red3", col.empirical = "navy", ...){
  g <- length(pi)
  if(n == floor(n) && n > 0 && length(alpha) == g && length(beta) == g &&
     length(gamma) == g && sum(pi) == 1 && min(c(pi, alpha, beta, gamma)) > 0 &&
     is.logical(plot.it) && is.logical(empirical)
  ){
    
    z = rmultinom(n = n, size = 1, pi)
    aux = rowSums(z)
    modal = dglindley_mix(moglindley_mix(pi, alpha, beta, gamma), pi, alpha, beta, gamma)
    
    sample = NULL
    for(j in 1:g){
      sample = c(sample, rglindley(aux[j], alpha[j], beta[j], gamma[j])$sample)
    }
    if(plot.it){
      d.breaks <- ceiling(nclass.Sturges(sample)*2.5)
      modal = max(modal, max(density(sample)$y))
      hist(sample, freq = F, border = "gray48",
           main = "Sampling distribution of X",xlab = "x",
           ylab = "Density",
           ylim = c(0, modal),
           if(any(names(list(...)) == "breaks") == FALSE){breaks = d.breaks}, ...)
      pop = function(x){dglindley_mix(x, pi, alpha, beta, gamma)}
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
      output = list(sample[,1], g, pi, alpha, beta, gamma, sample[,2], p)
      names(output) = c("sample", "g", "pi", "alpha", "beta", "gamma", "classification",
                        "plot")
    }
    else{
      output = list(sample[,1], g, pi, alpha, beta ,gamma, sample[,2])
      names(output) = c("sample", "g", "pi", "alpha", "beta", "gamma", "classification")
    }
    return(output)}
  else{
    stop("Error.")
  }
}
