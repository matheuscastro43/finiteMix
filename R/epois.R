epois = function(data, plot.it = TRUE, empirical = FALSE, 
                col.estimated = "orange", col.empirical = "navy", ...){
  if(is.numeric(data) && is.logical(plot.it) && is.logical(empirical)){
    data = sort(data)
    n = length(data)
    
    lambda = mean(data)
    
    lv = sum(dpois(data, lambda, log = TRUE))
    aic = 2 - 2*lv
    bic = log(n) - 2*lv
    
    if(plot.it == TRUE){
      d.breaks <- ceiling(nclass.Sturges(data)*2.5)
      modal = max(max(dpois(floor(lambda), lambda)),
                  hist(data, plot = FALSE, if(any(names(list(...)) == "breaks") == FALSE){
                    breaks = d.breaks}, ...)$density)
      hist(data, freq = F,border = "gray48",
           main = "Sampling distribution of X",xlab = "x",
           ylab = "Density",
           if(any(names(list(...)) == "breaks") == FALSE){breaks = 
             d.breaks}, ...)
      estimada = function(x){dpois(x, lambda)}
      x0 = min(data)
      x1 = max(data)
      par(new = T)
      plot(x0:x1, lapply(x0:x1, estimada), type = "h", col = col.estimated,
           lwd = 3, main = "",
           xlab = "", ylab = "", axes = F)
      if(empirical){
        lines(density(data),col = col.empirical,lwd = 3)
        legend("topright", legend=(c("Empirical", "Estimated")),
               fill=c(col.empirical, col.estimated),
               border = c(col.empirical, col.estimated), bty="n")
      }
      else{
        legend("topright", legend = "Estimated", fill = col.estimated,
               border = col.estimated, bty="n")
      }
      p <- recordPlot()
    }
    if(plot.it){
      output = list(lambda, lv, aic, bic, p)
      names(output) = c("lambda_hat", "logLik", "AIC", "BIC", "plot")}
    else{
      output = list(lambda, lv, aic, bic)
      names(output) = c("lambda_hat", "logLik", "AIC", "BIC")
    }
    return(output)
  }
}
