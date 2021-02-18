enorm = function(data, plot.it = TRUE, empirical = FALSE, 
                 col.estimated = "orange", col.empirical = "navy", ...){
  if(is.numeric(data) && is.logical(plot.it) && is.logical(empirical)){
    data = sort(data)
    n = length(data)
    
    medias = mean(data)
    dps = sd(data)
    
    lv = sum(log(dnorm(data, medias, dps)))
    
    if(plot.it == TRUE){
      modal = dnorm(medias, medias, dps)
      d.breaks = ceiling(nclass.Sturges(data)*2.5)
      modal = min(c(1, max(modal, hist(data, if(any(names(list(...)) ==
                                                    "breaks") == FALSE){
        breaks = d.breaks}, ...)$density)))
      hist(data, freq = F,border = "gray48",
           main = "Sampling distribution of X", xlab = "x",
           ylab = "Density",
           ylim = c(0, modal),
           if(any(names(list(...)) == "breaks") == FALSE){breaks = d.breaks},
           ...)
      
      estimada = function(x){dnorm(x, medias, dps)}
      curve(estimada, col = col.estimated, lwd = 3, add = T)
      if(empirical){
        lines(density(data),col = col.empirical, lwd = 3)
        legend("topright", legend=(c("Empirical", "Estimated")),
               fill=c(col.empirical, col.estimated), border = c(col.empirical,
                                                                col.estimated),
               bty="n")
      }
      else{
        legend("topright", legend = "Estimated", fill = col.estimated,
               border = col.estimated,
               bty="n")
      }
      p <- recordPlot()
    }
    if(plot.it){
      output = list(medias, dps, lv, p)
      names(output) = c("mu_hat", "sigma_hat", "logLik", "plot")}
    else{
      output = list(medias, dps, lv)
      names(output) = c("mu_hat", "sigma_hat", "logLik")
    }
    return(output)
  }
}
