eexp = function(data, plot.it = TRUE, empirical = FALSE, 
                 col.estimated = "orange", col.empirical = "navy", ...){
  if(is.numeric(data) && is.logical(plot.it) && is.logical(empirical)){
    data = sort(data)
    n = length(data)
    
    lambda = 1/mean(data)
    
    lv = sum(dexp(data, lambda, log = TRUE))
    
    if(plot.it){
      d.breaks = ceiling(nclass.Sturges(data)*2.5)
      modal = max(hist(data, plot = FALSE, if(any(names(list(...)) == "breaks") == FALSE){
        breaks = d.breaks}, ...)$density)
      hist(data, freq = F, border = "gray48",
           main = "Sampling distribution of X", xlab = "x",
           ylab = "Density",
           ylim = c(0, modal),
           if(any(names(list(...)) == "breaks") == FALSE){breaks = d.breaks},
           ...)
      
      estimada = function(x){dexp(x, lambda)}
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
      output = list(lambda, lv, p)
      names(output) = c("lambda_hat", "logLik", "plot")}
    else{
      output = list(lambda, lv)
      names(output) = c("lambda_hat", "logLik")
    }
    return(output)
  }
}
