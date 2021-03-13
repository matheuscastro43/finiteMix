egamma_mix = function(data, g, lim.em = 100, criteria = "dif.psi", plot.it = 
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
      est = egamma(data[k$cluster == j], plot.it = F)
      psi[1:3 > 1, j] = c(est$alpha_hat, est$beta_hat)
    }
    alphas = psi[2,]
    betas = psi[3,]
    medias = alphas * betas
    k = kmeans(data, g, centers = medias)
    psi[1, ] = pi = table(k$cluster)/n
    for(j in 1:g){
      est = egamma(data[k$cluster == j], plot.it = F)
      psi[1:3 > 1, j] = c(est$alpha_hat, est$beta_hat)
    }
    alphas = psi[2,]
    betas = psi[3,]
    
    count = 0
    L = function(i){dgamma_mix(data[i], pi, alphas, betas, log = TRUE)}
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
          Wij[i,j] <- as.numeric((pi[j]*dgamma(data[i], shape = alphas[j],
                                               scale = betas[j]))/
                                   sum((pi * dgamma(data[i], shape = alphas, 
                                                    scale = betas))))
        }
      }
      Wj <- colSums(Wij)
      pi <- 1/n * Wj
      pi = pi/sum(pi)
      
      Q = function(param){
        alphast = param[(1):(g)]
        betast = param[(g + 1):(2*g)]
        
        Qj = function(j){
          aux = 0
          for(i in 1:n){
            aux = aux + Wij[i,j] * log(pi[j] * dgamma(data[i], alphast[j],
                                                      scale = betast[j]))
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
          gr1 = gr1 + (Wij[i, ] * log(data[i]) - log(betast) * Wij[i, ] -
                         digamma(alphast) * Wij[i, ])
          gr2 = gr2 + (Wij[i, ] * data[i]/betast^2 - alphast/betast * Wij[i, ])
        }
        gr = c(gr1, gr2)
        return(-gr)
      }
      
      estim = optim(par = c(alphas, betas), fn = Q, method = "L-BFGS-B",
                    lower = rep(0.1, 3 * g),
                    upper = rep(Inf, 3 * g), gr = grQ, 
                    control = list(maxit = 10))
      alphas = estim$par[1:g]
      betas = estim$par[(g + 1):(2*g)]
      
      psi_new = matrix(c(pi, alphas, betas), 3, byrow = T)
      LF_new = sum(sapply(1:n, L))
      if(criteria == "dif.lh"){
        crit = LF_new - LF
        if((abs(crit) < 1*10^(-5))){cat("\n"); break}
        LF <- LF_new
      }
      else{
        crit = max(abs(psi - psi_new))
        if(any(is.na(crit))){
          alphas = start$alphas
          betas = start$betas
          medias = alphas * betas
          k = kmeans(data, g, centers = medias)
          pi = table(k$cluster)/n
          psi <- matrix(c(pi, alphas, betas), 3, byrow = T)
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
  p = prod(dim(matrix(0, 2, 2))) - 1
  aic = 2*p - 2*LF_new
  bic = p*log(n) - 2*LF_new
  if(plot.it == TRUE){
    d.breaks = ceiling(nclass.Sturges(data)*2.5)
    modal = 0
    for(i in 1:g){
      if(alphas[i] >= 1){
        modal[i] = max(dgamma(c((alphas[i]-1)*betas[i], (alphas[i])*
                                  betas[i]), alphas[i], scale = betas[i]))
      }else{
        U = modal[i] = 30
        while(modal[i] >= 0.9 * U){
          modal[i] = optimize(function(x) dgamma(x, alphas[i], 
                                                 scale = betas[i]),
                              interval = c(0, U), maximum = T)$maximum
          U = 2 * U
        }
        modal[i] = dgamma_mix(modal[i], pi, alphas, betas)
      }
      if(modal[i] > 10* dgamma_mix(1, pi, alphas, betas)){
        modal[i] = dgamma_mix(1, pi, alphas, betas)
      }
    }
    modal = min(c(1, max(modal, hist(data, plot = FALSE, 
                                     if(any(names(list(...)) == 
                                                  "breaks") == FALSE){breaks = 
                                                    d.breaks}, ...)$density)))
    hist(data, freq = F, border = "gray48",
         main = "Sampling distribution of X", xlab = "x",
         ylab = "Density",
         ylim = c(0, modal),
         if(any(names(list(...)) == "breaks") == FALSE){breaks = d.breaks}, ...)
    
    estimada = function(x){dgamma_mix(x, pi, alphas, betas)}
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
  medias = alphas * betas
  ordem = order(medias)
  medias = medias[ordem]
  pi = pi[ordem]
  alphas = alphas[ordem]
  betas = betas[ordem]
  gammal = function(x) {digamma(x) * gamma(x)}
  si = function(i){
    si = t(t(c((dgamma(data[i], alphas[-g], scale = betas[-g]) - 
                  dgamma(data[i], alphas[g], scale = betas[g]))/
                 dgamma_mix(data[i], pi, alphas, betas),
               pi * (exp(-data[i]/betas) * (betas^(-alphas) * 
                                              data[i]^(alphas - 1) * 
                                              (log(data[i]) - log(betas)) * 
                                              gamma(alphas) - betas^(-alphas) * 
                                              data[i]^(alphas - 1) * 
                                              gammal(alphas))/
                       (gamma(alphas))^2)/
                 dgamma_mix(data[i], pi, alphas, betas),
               pi * ((data[i]^(alphas - 1))/(gamma(alphas)) * 
                       ((-alphas) * betas^(-(alphas + 1)) * 
                          exp(-data[i]/betas) + betas^(-alphas) * 
                          exp(-data[i]/betas) * data[i]/betas^2))/
                 dgamma_mix(data[i], pi, alphas, betas))))
    rownames(si) = c(paste0("pi_", as.character(1:(g-1))), 
                     paste0("alpha_", as.character(1:(g))),
                     paste0("beta_", as.character(1:(g))))
    si %*% t(si)
  }
  se = sqrt(diag(solve(Reduce('+', sapply(1:n, si, simplify = FALSE)))))
  class = kmeans(data, centers = medias)$cluster
  if(plot.it){
    output = list(class, pi, alphas, betas, se, LF_new, 
                  aic, bic, count, p)
    names(output) = c("classification", "pi_hat", "alpha_hat", "beta_hat", 
                      "stde", "logLik", "AIC", "BIC", "EM-iterations", "plot")}
  else{
    output = list(class, pi, alphas, betas, se, LF_new, 
                  aic, bic, count)
    names(output) = c("classification", "pi_hat", "alpha_hat", "beta_hat",
                      "stde", "logLik", "AIC", "BIC", "EM-iterations")
  }
  return(output)
}
