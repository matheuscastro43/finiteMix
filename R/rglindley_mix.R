rglindley_mix = function(n, pi, alpha, beta, gamma, plot.it = TRUE, empirical = FALSE, 
                         col.pop = "red3", col.empirical = "navy", ...){
  g = length(pi)
  if(n == floor(n) && sum(pi) == 1 && min(c(pi, alpha, beta, gamma, n)) > 0 
     && length(alpha) == g && length(beta) == g && length(gamma) == g){
    
    z = rmultinom(n = n, size = 1, pi)
    aux = rowSums(z)
    modal = 0
    for(j in 1:g){
      if(alpha[j] >= 1){
        modal[j] = max(dglindley_mix(c((alpha[j]-1)*beta, (alpha[j])*beta), pi, alpha, beta, gamma))
      }else{
        U = modal[j] = 30
        while(modal[j] >= 0.9 * U){
          modal[j] = optimize(function(x) dglindley(x, alpha[j], beta[j], gamma[j]), interval = c(0, U), maximum = T)$maximum
          U = 2 * U
        }
        modal[j] = dglindley_mix(modal[j], pi, alpha, beta, gamma)
        if(modal[j] > 10* dglindley_mix(1, pi, alpha, beta, gamma)){
          modal[j] = dglindley_mix(1, pi, alpha, beta, gamma)
        }
      }
    }
    
    sample = NULL
    for(j in 1:g){
      sample = c(sample, rglindley(aux[j], alpha[j], beta[j], gamma[j], 
                                   plot.it = FALSE)$sample)
    }
    if(plot.it){
      d.breaks <- ceiling(nclass.Sturges(sample)*2.5)
      modal = min(c(1, max(modal, hist(sample, plot = FALSE,
                                       if(any(names(list(...)) == "breaks") ==
                                          FALSE){
                                         breaks = d.breaks}, ...)$density)))
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
  else stop("The parametric space must be respected.")
}
