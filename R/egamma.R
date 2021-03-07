egamma = function(data, plot.it = TRUE, empirical = FALSE, 
                     col.estimated = "orange", col.empirical = "navy", ...){
  if((is.numeric(data) || is.numeric(data$sample)) && is.logical(plot.it) &&
     is.logical(empirical)){
    if(is.list(data)){data = data$sample}
    data = sort(data)
    n = length(data)
    
    alphas = mean(data)^2/var(data)
    betas = var(data)/mean(data)
    
    LV = function(Psi, x){
      alpha = Psi[1]
      beta = Psi[2]
      lv = (alpha - 1) * sum(log(x)) - sum(x)/beta - n * alpha * log(beta) -
        n * log(gamma(alpha))
      if(lv == -Inf) return(.Machine$double.xmax/1e+08)
      if(lv == Inf) return(-.Machine$double.xmax/1e+08)
      return(-lv)
    }
    
    grr = function(Psi, x){
      alpha = Psi[1]
      beta = Psi[2]
      -c(sum(log(x)) - n * digamma(alpha) - n * log(beta),
         sum(x)/beta^2 - n * alpha/beta)
    }
    b = optim(par = c(alphas, betas), fn = LV, lower = c(1e-04, 1e-04), 
              upper = c(Inf, Inf), method = "L-BFGS-B", x = data, gr = grr)
    
    alpha = b$par[1]
    beta = b$par[2]
    LF = -b$value
    
    if(alpha >= 1){
      modal = max(dgamma(c((alpha-1)*beta, (alpha)*beta), alpha, scale = beta))
    }else{
      U = modal = 30
      while(modal >= 0.9 * U){
        modal = optimize(function(x) dgamma(x, alpha, scale = beta),
                         interval = c(0, U), maximum = T)$maximum
        U = 2 * U
      }
      modal = dgamma(modal, alpha, scale = beta)
      if(modal > 10* dgamma(1, alpha, scale = beta)){
        modal = dgamma(1, alpha, scale = beta)
      }
    }
    
    if(plot.it == TRUE){
      d.breaks = ceiling(nclass.Sturges(data)*2.5)
      modal = min(c(1, max(modal, hist(data, plot = FALSE, if(any(names(list(...)) == "breaks") == FALSE){
        breaks = d.breaks}, ...)$density)))
      hist(data,freq = F,border = "gray48",
           main = "Sampling distribution of X", xlab = "x",
           ylab = "Density",
           ylim = c(0, modal),
           if(any(names(list(...)) == "breaks") == FALSE){breaks = d.breaks}, ...)
      
      estimada = function(x){dgamma(x, alpha, scale = beta)}
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
