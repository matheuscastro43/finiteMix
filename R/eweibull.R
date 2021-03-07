eweibull = function(data, plot.it = TRUE, empirical = FALSE, 
                    col.estimated = "orange", col.empirical = "navy", ...){
  if((is.numeric(data) || is.numeric(data$sample)) && is.logical(plot.it) &&
     is.logical(empirical)){
    if(is.list(data)){data = data$sample}
    data = sort(data)
    n = length(data)
    
    medias = mean(data)
    vars = var(data)
    
    shapeMoments <- function(shape){
      log(gamma(1 + 2/shape)) - 2*log(gamma(1 + 1/shape)) - log(vars + medias^2) +
        2*log(medias)
    }
    shapes = uniroot(shapeMoments, c(0.1, 100))$root
    scales = medias/gamma(1 + 1/shapes)
    
    LV <- function(Psi, x){
      shape = Psi[1]
      scale = Psi[2]
      lv = n*log(shape) - n*log(scale) + (shape - 1)*sum((log(x))) -
        n*(shape - 1)*log(scale) - sum((x/scale)^shape)
      if(lv == -Inf) return(.Machine$double.xmax/1e+08)
      if(lv == Inf) return(-.Machine$double.xmax/1e+08)
      return(-lv)
    }
    
    grr = function(Psi, x){
      shape = Psi[1]
      scale = Psi[2]
      -c(n/shape + sum(log(x)) - n*log(scale) - 
           sum(log(x/scale) * (x/scale)^shape), 
         -n/scale - n*(shape - 1)/scale + sum(shape/scale * (x/scale)^shape))
    }
    b = optim(par = c(1, 9), fn = LV, lower = c(1e-04, 1e-04), 
              upper = c(Inf, Inf), method = "L-BFGS-B", x = data, gr = grr)
    
    alpha = b$par[1]
    beta = b$par[2]
    LF = -b$value
    
    if(alpha >= 1){
      modal = dweibull(beta * ((alpha - 1)/alpha)^(1/beta), alpha, beta)
    }else{
      modal = dweibull(.1, alpha, beta)
    }
    
    if(plot.it == TRUE){
      d.breaks = ceiling(nclass.Sturges(data)*2.5)
      modal = max(modal, hist(data, plot = FALSE, if(any(names(list(...)) == "breaks") == FALSE){
        breaks = d.breaks}, ...)$density)
      hist(data,freq = F,border = "gray48",
           main = "Sampling distribution of X", xlab = "x",
           ylab = "Density",
           ylim = c(0, modal),
           if(any(names(list(...)) == "breaks") == FALSE){breaks = d.breaks}, ...)
      
      estimada = function(x){dweibull(x, alpha, beta)}
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
    ordem = order(alpha)
    if(plot.it){
      output = list(alpha[ordem], beta[ordem], LF, p)
      names(output) = c("alpha_hat", "beta_hat", "logLik", "plot")}
    else{
      output = list(alpha[ordem], beta[ordem], LF)
      names(output) = c("alpha_hat", "beta_hat", "logLik")
    }
    return(output)
  }
}
