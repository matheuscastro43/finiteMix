eglindley = function(data, plot.it = TRUE, empirical = FALSE, 
                     col.estimated = "orange", col.empirical = "navy", ...){
  if((is.numeric(data) || is.numeric(data$sample)) && is.logical(plot.it) &&
     is.logical(empirical)){
    if(is.list(data)){data = data$sample}
    data <- sort(data)
    n <- length(data)
    
    invisible(utils::capture.output(b <- bmixture::bmixgamma(data, k = 2)))
    alphas = mean(apply(b$alpha_sample, 2, mean)) - .5
    betas = mean(apply(b$beta_sample, 2, mean))
    gammas = (b$pi_sample[1]^(-1) -1)/betas
    
    LV = function(Psi, x){
      alpha = Psi[1]
      beta = Psi[2]
      gamma = Psi[3]
      lv = (alpha - 1) * sum(log(x)) + sum(log(alpha + gamma*x)) - sum(x)/beta -
        n * (alpha * log(beta) + log(beta*gamma + 1) + log(gamma(alpha + 1)))
      return(-lv)
    }
    
    grr = function(Psi, x){
      alpha = Psi[1]
      beta = Psi[2]
      gamma = Psi[3]
      -c(sum(log(x)) + sum(1/(alpha + gamma * x)) - n * log(beta) - n * 
                                 digamma(alpha + 1),
             sum(x)/beta^2 - n*alpha/beta - n * gamma/(beta * gamma + 1),
             sum(x/(alpha + gamma * x)) - n * beta/(beta*gamma + 1))
    }
    b = optim(par = c(alphas, betas, gammas), fn = LV, lower = c(1e-04, 1e-04, 1e-04), 
              upper = c(Inf, Inf, Inf), method = "L-BFGS-B", x = data, gr = grr)
    
    alpha = b$par[1]
    beta = b$par[2]
    gamma = b$par[3]
    
    if(alpha >= 1){
      modal = max(dglindley(c((alpha-1)*beta, (alpha)*beta), alpha, beta,
                            gamma))
    }else{
      U = modal = 30
      while(modal >= 0.9 * U){
        modal = optimize(function(x) dglindley(x, alpha, beta, gamma),
                         interval = c(0, U), maximum = T)$maximum
        U = 2 * U
      }
      modal = dglindley(modal, alpha, beta, gamma)
      if(modal > 10* dglindley(1, alpha, beta, gamma)){
        modal = dglindley(1, alpha, beta, gamma)
      }
    }
    
    if(plot.it == TRUE){
      d.breaks = ceiling(nclass.Sturges(data)*2.5)
      modal = min(c(1, max(modal, hist(data, if(any(names(list(...)) == "breaks") == FALSE){
        breaks = d.breaks}, ...)$density)))
      hist(data,freq = F,border = "gray48",
           main = "Sampling distribution of X", xlab = "x",
           ylab = "Density",
           ylim = c(0, modal),
           if(any(names(list(...)) == "breaks") == FALSE){breaks = d.breaks}, ...)
      
      estimada = function(x){dglindley(x, alpha, beta, gamma)}
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
      output = list(alpha[ordem], beta[ordem], gamma[ordem], p)
      names(output) = c("alpha_hat", "beta_hat", "gamma_hat", "plot")}
    else{
      output = list(alpha[ordem], beta[ordem], gamma[ordem])
      names(output) = c("alpha_hat", "beta_hat", "gamma_hat")
    }
    return(output)
  }
}
