eweibull_mix = function(data, g, lim.em = 100, criteria = "dif.psi", 
                        epsilon = 1e-05, plot.it = TRUE, empirical = FALSE, 
                        col.estimated = "orange", col.empirical = "navy", ...){
  if((is.numeric(data) || is.numeric(data$sample)) && g == floor(g) && g > 1 &&
     is.logical(plot.it) && is.logical(empirical) &&
     (criteria == "dif.lh" || criteria == "dif.psi")){
    if(is.list(data)){data = data$sample}
    data = sort(data)
    n = length(data)
    
    psi = matrix(0, ncol = g, nrow = 3)
    k = kmeans(data, g)
    for(j in 1:g){
      est = eweibull(data[k$cluster == j], plot.it = F)
      psi[1:3 > 1, j] = c(est$alpha_hat, est$beta_hat)
    }
    alphas = psi[2,]
    betas = psi[3,]
    medias = betas * gamma(1 + 1/alphas)
    k = kmeans(data, g, centers = medias)
    psi[1, ] = pi = table(k$cluster)/n
    for(j in 1:g){
      est = eweibull(data[k$cluster == j], plot.it = F)
      psi[1:3 > 1, j] = c(est$alpha_hat, est$beta_hat)
    }
    alphas = psi[2,]
    betas = psi[3,]
    
    count = 0
    L = function(i){dweibull_mix(data[i], pi, alphas, betas, log = TRUE)}
    LF = sum(sapply(1:n, L))
    while(T){
      progress <- function (x, max = lim.em) {
        percent <- x / max * 100
        cat(sprintf('\r[%-50s] %d%%',
                    paste(rep('=', percent / 2), collapse = ''),
                    floor(percent)))
        if (x == max)
          cat('\n')
      }
      if(count == 0)
        cat("Limit of EM Iterations (", lim.em ,"): \n", sep = "")
      progress(count)
      Wij = matrix(0, nrow = n, ncol = g)
      for(i in 1:n){
        for(j in 1:g){
          Wij[i,j] <- as.numeric((pi[j]*dweibull(data[i], alphas[j], 
                                                 betas[j]))/
                                   sum((pi * dweibull(data[i], alphas, 
                                                      betas))))
        }
      }
      Wj <- colSums(Wij)
      pi <- 1/n * Wj
      pi = pi/sum(pi)
      if(any(is.nan(pi))){
        medias = betas * gamma(1 + 1/alphas)
        psi[1, ] = pi = table(kmeans(data, g, centers = medias)$cluster)/n
      }
      
      Q = function(param){
        alphast = param[(1):(g)]
        betast = param[(g + 1):(2*g)]
        
        
        aux = 0
        for(i in 1:n){
          aux2 = Wij[i,j] * dweibull_mix(data[i], pi, alphast, betast, 
                                         log = TRUE)
          if(is.nan(aux2)){next}
          if(aux2 == Inf){aux = Inf; break
          }else{
            if(aux2 == -Inf){aux = -Inf; break}
            else{
              aux = aux + aux2}}
        }
        
        if(aux == -Inf) return(.Machine$double.xmax/1e+100)
        if(aux == Inf) return(-.Machine$double.xmax/1e+100)
        return(-aux)
      }
      grQ = function(param){
        alphast = param[(1):(g)]
        betast = param[(g + 1):(2*g)]
        
        gr1 = gr2 = rep(0, g)
        for(i in 1:n){
          gr1 = gr1 + (Wij[i, ]/alphast + Wij[i, ]*log(data[i]) -
                         Wij[i, ]*log(betast) - Wij[i, ] * 
                         log(data[i]/betast*(data[i]/betast)^alphast))
          gr2 = gr2 + (-Wij[i, ]/betast - Wij[i, ]*(alphast - 1)/betast +
                         alphast/betast * Wij[i, ]*(data[i]/betast)^alphast)
        }
        gr = c(gr1, gr2)
        return(-gr)
      }
      estim = NULL
      while(is.null(estim)){
        estim = tryCatch(optim(par = c(alphas, betas), fn = Q, 
                               method = "L-BFGS-B", lower = rep(1e-01, 2 * g),
                               upper = rep(Inf, 2 * g), gr = grQ),
                         error = function(e){NULL})
        
        if(is.null(estim)){
          alphas = alphas + 0.1
          betas = betas + 0.1
          count = max(c(0, count - 1))
        }
      }
      alphas = estim$par[1:g]
      betas = estim$par[(g + 1):(2*g)]
      
      psi_new = matrix(c(pi, alphas, betas), 3, byrow = T)
      LF_new = sum(sapply(1:n, L))
      if(criteria == "dif.lh"){
        crit = LF_new - LF
        if((abs(crit) < epsilon)){cat("\n"); break}
        LF <- LF_new
      }
      else{
        crit = max(abs(psi - psi_new))
        if(any(is.na(crit))){
          k = kmeans(data, g)
          for(j in 1:g){
            est = eweibull(data[k$cluster == j], plot.it = F)
            psi[1:3 > 1, j] = c(est$alpha_hat, est$beta_hat)
          }
          alphas = psi[2,]
          betas = psi[3,]
          medias = betas * gamma(1 + 1/alphas)
          k = kmeans(data, g, centers = medias)
          psi[1, ] = pi = table(k$cluster)/n
          for(j in 1:g){
            est = eweibull(data[k$cluster == j], plot.it = F)
            psi[1:3 > 1, j] = c(est$alpha_hat, est$beta_hat)
          }
          alphas = psi[2,]
          betas = psi[3,]
          next
        }
        if(crit < epsilon) {cat("\n"); break}
        psi = psi_new
      }
      count = count + 1
      if(count >= lim.em){
        progress(count)
        message("\nLimit of Iterations reached!")
        break
      }
    }
  }
  p = 3*g - 1
  aic = 2*p - 2*LF_new
  bic = p*log(n) - 2*LF_new
  if(plot.it == TRUE){
    d.breaks = ceiling(nclass.Sturges(data)*2.5)
    modal = max(dweibull_mix(c(0.1, 
                               betas * ((alphas - 1)/alphas)^(alphas - 1)),
                             pi, alphas, betas))
    
    
    modal = max(modal, hist(data, plot = FALSE, if(any(names(list(...)) == 
                                                       "breaks") == FALSE){
      breaks = d.breaks}, ...)$density)
    hist(data, freq = F, border = "gray48",
         main = "Sampling distribution of X", xlab = "x",
         ylab = "Density",
         ylim = c(0, modal),
         if(any(names(list(...)) == "breaks") == FALSE){breaks = d.breaks}, ...)
    
    estimada = function(x){dweibull_mix(x, pi, alphas, betas)}
    curve(estimada, col = col.estimated, lwd = 3, add = T)
    if(empirical){
      lines(density(data),col = col.empirical, lwd = 3)
      legend("topright", legend=(c("Empirical", "Estimated")),
             fill=c(col.empirical, col.estimated), border = c(col.empirical, 
                                                              col.estimated),
             bty="n")
    }
    else{
      legend("topright", legend = "Estimated", fill = col.estimated, border = 
               col.estimated, bty="n")
    }
    p <- recordPlot()
  }
  medias = betas * gamma(1 + 1/alphas)
  ordem = order(medias)
  medias = medias[ordem]
  pi = pi[ordem]
  alphas = alphas[ordem]
  betas = betas[ordem]
  si = function(i){
    si = t(t(c((dweibull(data[i], alphas[-g], betas[-g]) - 
                  dweibull(data[i], alphas[g], betas[g]))/
                 dweibull_mix(data[i], pi, alphas, betas),
               pi * ((1/betas * (data[i]/betas)^(alphas - 1) + alphas/betas *
                        (data[i]/betas)^(alphas - 1) * log(data[i]/betas)) * 
                       exp(-(data[i]/betas)^alphas) - alphas/betas * 
                       (data[i]/betas)^(alphas - 1) * 
                       exp(-(data[i]/betas)^(alphas)) * 
                       (data[i]/betas)^alphas * log(data[i]/betas))/
                 dweibull_mix(data[i], pi, alphas, betas),
               pi * ((alphas/betas * data[i]^(alphas - 1) * (1 - alphas) * 
                        betas^(-alphas) - alphas/betas^2 * 
                        (data[i]/betas)^(alphas - 1)) * 
                       exp(-(data[i]/betas)^(alphas)) + alphas/betas * 
                       (data[i]/betas)^(alphas - 1) * 
                       exp(-(data[i]/betas)^(alphas)) * alphas * 
                       data[i]^(alphas) * betas^(-(alphas + 1)))/
                 dweibull_mix(data[i], pi, alphas, betas))))
    rownames(si) = c(paste0("pi_", as.character(1:(g-1))),
                     paste0("alpha_", as.character(1:g)),
                     paste0("beta_", as.character(1:g)))
    si %*% t(si)
  }
  se = sqrt(diag(solve(Reduce('+', sapply(1:n, si, simplify = FALSE)))))
  class = kmeans(data, centers = medias)$cluster
  if(plot.it){
    output = list(class, pi, alphas, betas, 
                  se, LF_new, aic, bic, count, p)
    names(output) = c("classification", "pi_hat", "alpha_hat", "beta_hat", 
                      "stde", "logLik", "AIC", "BIC", "EM_iterations", "plot")}
  else{
    output = list(class, pi, alphas, betas,
                  se, LF_new, aic, bic, count)
    names(output) = c("classification", "pi_hat", "alpha_hat", "beta_hat",
                      "stde", "logLik", "AIC", "BIC", "EM_iterations")
  }
  return(output)
}
