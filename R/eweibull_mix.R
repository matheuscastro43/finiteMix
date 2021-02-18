eweibull_mix = function(data, g, lim.em = 100, criteria = "dif.psi", plot.it = 
                          TRUE, empirical = FALSE, col.estimated = "orange", 
                        col.empirical = "navy", ...){
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
    
    L = function(I){sum(log(pi[I] * dweibull(data[which(k$cluster == I)], 
                                             alphas[I], betas[I])))}
    LF = sum(as.numeric(lapply(1:g, L))); count = 0
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
        cat("Limit of EM Interactions (", lim.em ,"): \n", sep = "")
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
      
      Q = function(param){
        alphast = param[(1):(g)]
        betast = param[(g + 1):(2*g)]
        
        Qj = function(j){
          aux = 0
          for(i in 1:n){
            aux = aux + Wij[i,j] * (log(pi[j]) + log(dweibull(data[i], alphast[j],
                                                        betast[j])))
            if(aux == -Inf) return(aux = -.Machine$double.xmax/1e+08)
            if(aux == Inf) return(aux = .Machine$double.xmax/1e+08)
          }
          return(aux)
        }
        q = sum(sapply(1:g, Qj))
        if(q == -Inf) return(.Machine$double.xmax/1e+08)
        if(q == Inf) return(-.Machine$double.xmax/1e+08)
        return(-q)
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
      
      estim = optim(par = c(alphas, betas), fn = Q, method = "L-BFGS-B",
                    lower = rep(1e-04, 2 * g),
                    upper = rep(Inf, 2 * g), gr = grQ)
      alphas = estim$par[1:g]
      betas = estim$par[(g + 1):(2*g)]
      
      psi_new = matrix(c(pi, alphas, betas), 3, byrow = T)
      LF_new = sum(as.numeric(lapply(1:g, L)))
      if(criteria == "dif.lh"){
        crit = LF_new - LF
        if((abs(crit) < 1*10^(-5))){cat("\n"); break}
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
        if(crit < 1*10^(-5)) {cat("\n"); break}
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
  if(plot.it == TRUE){
    d.breaks = ceiling(nclass.Sturges(data)*2.5)
    modal = max(dweibull_mix(c(0.1, 
                               betas * ((alphas - 1)/alphas)^(alphas - 1)),
                             pi, alphas, betas))
    
    
    modal = max(modal, hist(data, if(any(names(list(...)) == 
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
  class = kmeans(data, centers = medias[ordem])$cluster
  if(plot.it){
    output = list(class, pi[ordem], alphas[ordem], betas[ordem], 
                  LF_new, count, p)
    names(output) = c("classification", "pi_hat", "alpha_hat", "beta_hat", 
                      "logLik", "EM-interactions", "plot")}
  else{
    output = list(class, pi[ordem], alphas[ordem], betas[ordem],
                  LF_new, count)
    names(output) = c("classification", "pi_hat", "alpha_hat", "beta_hat",
                      "logLik", "EM-interactions")
  }
  return(output)
}
