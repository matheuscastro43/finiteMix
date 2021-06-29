eglindley_mix = function(data, g, lim.em = 100, criteria = "dif.psi", 
                         epsilon = 1e-05, plot.it = TRUE, empirical = FALSE, 
                         col.estimated = "orange", col.empirical = "navy", ...){
  if((is.numeric(data) || is.numeric(data$sample)) && g == floor(g) && g > 1 &&
     is.logical(plot.it) && is.logical(empirical) &&
     (criteria == "dif.lh" || criteria == "dif.psi")){
    if(is.list(data)){data = data$sample}
    data = sort(data)
    n = length(data)
    
    psi = matrix(0, ncol = g, nrow = 4)
    k = kmeans(data, g)
    for(j in 1:g){
      est = eglindley(data[k$cluster == j], plot.it = F)
      psi[1:4 > 1, j] = c(est$alpha_hat, est$beta_hat, est$gamma_hat)
    }
    alphas = psi[2,]
    betas = psi[3,]
    gammas = psi[4,]
    medias = alphas * betas + (betas^2 * gammas)/(1 + betas * gammas)
    k = kmeans(data, g, centers = medias)
    psi[1, ] = pi = table(k$cluster)/n
    for(j in 1:g){
      est = eglindley(data[k$cluster == j], plot.it = F)
      psi[1:4 > 1, j] = c(est$alpha_hat, est$beta_hat, est$gamma_hat)
    }
    alphas = psi[2,]
    betas = psi[3,]
    gammas = psi[4,]

    L = function(i){dglindley_mix(data[i], pi, alphas, betas, gammas, log = TRUE)}
    LF = sum(sapply(1:n, L)); count = 0
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
          Wij[i,j] <- as.numeric((pi[j]*dglindley(data[i], alpha = alphas[j],
                                                  beta = betas[j], 
                                                  gamma = gammas[j]))/
                                   sum((pi * dglindley(data[i], alpha = alphas, 
                                                       beta = betas,
                                                       gamma = gammas))))
        }
      }
      Wj <- colSums(Wij)
      pi <- 1/n * Wj
      pi = pi/sum(pi)
      
      Q = function(param){
        alphast = param[(1):(g)]
        betast = param[(g + 1):(2*g)]
        gammast = param[(2*g + 1):(3*g)]
        
        Qj = function(j){
          aux = 0
          for(i in 1:n){
            aux = aux + Wij[i,j] * log(pi[j] * dglindley(data[i], alphast[j],
                                                         betast[j], gammast[j]))
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
        gammast = param[(2*g + 1):(3*g)]
        
        gr1 = gr2 = gr3 = rep(0, g)
        for(i in 1:n){
          gr1 = gr1 + (Wij[i,] * log(data[i]) + Wij[i,]/
                         (alphast + gammast * data[i]) - log(betast) * 
                         Wij[i,] - digamma(alphast + 1) * Wij[i, ])
          gr2 = gr2 + ((Wij[i, ] * data[i])/betast^2 - alphast/betast * 
                         Wij[i,] - gammast/(betast * gammast + 1) * Wij[i, ])
          gr3 = gr3 + ((Wij[i, ] * data[i])/(alphast + gammast * data[i]) - 
                         betast/(betast * gammast + 1) * Wij[i,])
        }
        gr = c(gr1, gr2, gr3)
        return(-gr)
      }
      
      estim = optim(par = c(alphas, betas, gammas), fn = Q, method = "L-BFGS-B",
                    lower = rep(1e-04, 3 * g),
                    upper = rep(Inf, 3 * g), gr = grQ)
      alphas = estim$par[1:g]
      betas = estim$par[(g + 1):(2*g)]
      gammas = estim$par[(2*g + 1):(3*g)]
      
      psi_new = matrix(c(pi, alphas, betas, gammas), 4, byrow = T)
      LF_new = sum(sapply(1:n, L))
      if(criteria == "dif.lh"){
        crit = LF_new - LF
        if((abs(crit) < epsilon)){cat("\n"); break}
        LF <- LF_new
      }
      else{
        crit = max(abs(psi - psi_new))
        if(any(is.na(crit))){
          alphas = start$alphas
          betas = start$betas
          gammas = start$gammas
          medias = alphas * betas + (betas^2 * gammas)/(1 + betas * gammas)
          k = kmeans(data, g, centers = medias)
          pi = table(k$cluster)/n
          psi <- matrix(c(pi, alphas, betas, gammas), 4, byrow = T)
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
  p = 4*g - 1
  aic = 2*p - 2*LF_new
  bic = p*log(n) - 2*LF_new
  if(plot.it == TRUE){
    d.breaks = ceiling(nclass.Sturges(data)*2.5)
    modal = 0
    for(i in 1:g){
      if(alphas[i] >= 1){
        modal[i] = max(dglindley(c((alphas[i]-1)*betas[i], (alphas[i])*
                                     betas[i]), alphas[i], betas[i], gammas[i]))
      }else{
        U = modal[i] = 30
        while(modal[i] >= 0.9 * U){
          modal[i] = optimize(function(x) dglindley(x, alphas[i], betas[i],
                                                    gammas[i]),
                              interval = c(0, U), maximum = T)$maximum
          U = 2 * U
        }
      }
      modal[i] = dglindley_mix(modal[i], pi, alphas, betas, gammas)
      if(modal[i] > 10* dglindley_mix(1, pi, alphas, betas, gammas)){
        modal[i] = dglindley_mix(1, pi, alphas, betas, gammas)
      }
    }
    modal = min(c(1, max(modal, hist(data, plot = FALSE, if(any(names(list(...)) == 
                                                  "breaks") == FALSE){
      breaks = d.breaks}, ...)$density)))
    hist(data, freq = F, border = "gray48",
         main = "Sampling distribution of X", xlab = "x",
         ylab = "Density",
         ylim = c(0, modal),
         if(any(names(list(...)) == "breaks") == FALSE){breaks = d.breaks}, ...)
    
    estimada = function(x){dglindley_mix(x, pi, alphas, betas, gammas)}
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
  medias = alphas * betas + (betas^2 * gammas)/(1 + betas * gammas)
  ordem = order(medias)
  medias = medias[ordem]
  pi = pi[ordem]
  alphas = alphas[ordem]
  betas = betas[ordem]
  gammas = gammas[ordem]
  gammal = function(x) {digamma(x) * gamma(x)}
  si = function(i){
    si = t(t(c((dglindley(data[i], alphas[-g], betas[-g], gammas[-g]) - 
                  dglindley(data[i], alphas[g], betas[g], gammas[g]))/
                 dglindley_mix(data[i], pi, alphas, betas, gammas),
               pi * ((exp(-data[i]/betas)/(betas * gammas + 1)) * 
                       (((data[i]^(alphas - 1) * log(data[i]) * 
                            (alphas + gammas * 
                               data[i]) + data[i]^(alphas - 1)) * 
                           betas^(alphas) * gamma(alphas + 1) - 
                           data[i]^(alphas - 1) * (alphas + gammas * data[i]) * 
                           (betas^alphas * log(betas) * gamma(alphas + 1) + 
                              betas^alphas * gammal(alphas + 1)))/
                          (betas^alphas * gamma(alphas + 1))^2))/
                 dglindley_mix(data[i], pi, alphas, betas, gammas),
               pi * ((data[i]^(alphas - 1) * (alphas + gammas * data[i]))/
                       (gamma(alphas + 1)) * (exp(-data[i]/betas) * 
                                                (data[i]/betas^2) * 
                                                betas^alphas * 
                                                (betas * gammas + 1) - 
                                                exp(-data[i]/betas) * 
                                                (gammas*(alphas + 1) * 
                                                   betas^(alphas) + alphas * 
                                                   betas^(alphas - 1)))/
                       (betas^alphas * (betas * gammas + 1))^2)/
                 dglindley_mix(data[i], pi, alphas, betas, gammas),
               pi * ((data[i]^(alphas - 1) * exp(-data[i]/betas))/
                       (betas^alphas * gamma(alphas + 1)) * 
                       (data[i] - alphas * betas)/(betas * gammas + 1)^2)/
                 dglindley_mix(data[i], pi, alphas, betas, gammas))))
    rownames(si) = c(paste0("pi_", as.character(1:(g-1))),
                     paste0("alpha_", as.character(1:(g))),
                     paste0("beta_", as.character(1:(g))),
                     paste0("gamma_", as.character(1:(g))))
    si %*% t(si)
  }
  #se = sqrt(diag(solve(Reduce('+', sapply(1:n, si, simplify = FALSE)))))
  class = kmeans(data, centers = medias)$cluster
  if(plot.it){
    output = list(class, pi, alphas, betas, gammas, LF_new, aic, bic, count, p)
    names(output) = c("classification", "pi_hat", "alpha_hat", "beta_hat", 
                      "gamma_hat", "logLik", "AIC", "BIC", 
                      "EM_iterations", "plot")}
  else{
    output = list(class, pi, alphas, betas, gammas, LF_new, aic, bic, count)
    names(output) = c("classification", "pi_hat", "alpha_hat", "beta_hat",
                      "gamma_hat", "logLik", "AIC", "BIC", 
                      "EM_iterations")
  }
  return(output)
}
