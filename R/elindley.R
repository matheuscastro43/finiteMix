elindley = function(data, plot.it = TRUE, empirical = FALSE, 
                  col.estimated = "orange", col.empirical = "navy", ...){
  if((is.numeric(data) || is.numeric(data$sample)) && is.logical(plot.it) &&
     is.logical(empirical)){
    if(is.list(data)){data = data$sample}
    data = sort(data)
    n = length(data)
    
    betas = var(data)/mean(data)
    
    LV = function(Psi, x){
      beta = Psi
      lv = sum(log(x + 1)) - sum(x)/beta - n * log(beta) -
        n * log(beta + 1)
      if(lv == -Inf) return(.Machine$double.xmax/1e+08)
      if(lv == Inf) return(-.Machine$double.xmax/1e+08)
      return(-lv)
    }
    
    grr = function(Psi, x){
      beta = Psi
      -(-n/beta - n/(beta + 1) + sum(x)/beta^2)
    }
    b = optim(par = betas, fn = LV, lower = 1e-04, upper = Inf,
              method = "L-BFGS-B", x = data, gr = grr)
    
    beta = b$par
    
      modal = max(dlindley(c(0, beta), beta))
    
    if(plot.it == TRUE){
      d.breaks = ceiling(nclass.Sturges(data)*2.5)
      modal = min(c(1, max(modal, hist(data, if(any(names(list(...)) ==
                                                    "breaks") == FALSE){
        breaks = d.breaks}, ...)$density)))
      hist(data,freq = F,border = "gray48",
           main = "Sampling distribution of X", xlab = "x",
           ylab = "Density",
           ylim = c(0, modal),
           if(any(names(list(...)) == "breaks") == FALSE){breaks = d.breaks}, ...)
      
      estimada = function(x){dlindley(x, beta)}
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
    ordem = order(beta)
    if(plot.it){
      output = list(beta[ordem], p)
      names(output) = c("beta_hat", "plot")}
    else{
      output = list(beta[ordem])
      names(output) = c("beta_hat")
    }
    return(output)
  }
}
