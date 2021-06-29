enorm_mix = function(data, g, lim.em = 100, criteria = "dif.psi",
                     epsilon = 1e-05,plot.it = TRUE, empirical = FALSE, 
                     col.estimated = "orange", col.empirical = "navy", ...){
  if((is.numeric(data) || is.numeric(data$sample)) && g == floor(g) && g > 1 &&
     is.logical(plot.it) && is.logical(empirical) &&
     (criteria == "dif.lh" || criteria == "dif.psi")){
    if(is.list(data)){data = data$sample}
    data <- sort(data)
    n <- length(data)
    
    k <- kmeans(data, g)
    pi <- table(k$cluster)/n
    medias <- tapply(data, k$cluster, mean)
    dps <- tapply(data, k$cluster, sd)
    psi <- matrix(c(pi, medias, dps), 3, byrow = T)
    
    count = 0
    L <- function(i){dnorm_mix(data[i], pi, medias, dps, log = TRUE)}
    LF <- sum(sapply(1:n, L))
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
      Wij <- matrix(0, nrow = n, ncol = g)
      for(i in 1:n){
        for(j in 1:g){
          Wij[i,j] <- as.numeric((pi[j]*dnorm(data[i], mean = medias[j], 
                                              sd = dps[j]))/
                                   sum((pi * dnorm(data[i], mean = medias, 
                                                   sd = dps))))
        }
      }
      Wj <- colSums(Wij)
      pi <- 1/n * Wj
      pi = pi/sum(pi)
      medias <- as.numeric(lapply(1:g, function(j){
        aux = 0; for(i in 1:n){aux = aux + (data[i]*Wij[i,j])/(Wj[j])};
      return(aux)}))
      dps <- sqrt(as.numeric(lapply(1:g, function(j){
        aux = 0; for(i in 1:n){
          aux = aux + ((data[i]-medias[j])^2*Wij[i,j])/(Wj[j])};
      return(aux)})))
      psi_new <- matrix(c(pi, medias, dps), 3, byrow = T)
      LF_new <- sum(sapply(1:n, L))
      if(criteria == "dif.lh"){
        crit <- LF_new - LF
        if((abs(crit) < epsilon)){cat("\n"); break}
        LF <- LF_new
      }
      else{
        crit = max(abs(psi - psi_new))
        if(any(is.na(crit))){
          k <- kmeans(data, g)
          pi <- table(k$cluster)/n
          medias <- tapply(data, k$cluster, mean)
          dps <- tapply(data, k$cluster, sd)
          psi <- matrix(c(pi, medias, dps), 3, byrow = T)
          next
        }
        if(crit < epsilon){cat("\n"); break}
        psi <- psi_new
      }
      count = count + 1
      if(count >= lim.em){
        progress(count)
        message("\nLimit of Iterations reached!")
        break
      }
    }
    p = 3*g - 1
    aic = 2*p - 2*LF_new
    bic = p*log(n) - 2*LF_new
    if(plot.it){
      d.breaks <- ceiling(nclass.Sturges(data)*2.5)
      modal <- max(dnorm_mix(medias, pi, medias, dps))
      hist(data,freq = F,border = "gray48",
           main = "Sampling distribution of X", xlab = "x",
           ylab = "Density",
           ylim = c(0, modal),
           if(any(names(list(...)) == "breaks") == FALSE){breaks = 
             d.breaks}, ...)
      
      estimada = function(x){dnorm_mix(x, pi, medias, dps)}
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
    si = function(i){
      ordem = order(medias)
      pi = pi[ordem]
      medias = medias[ordem]
      dps = dps[ordem]
      si = t(t(c((dnorm(data[i], medias[-g], dps[-g]) - 
                    dnorm(data[i],medias[g], dps[g]))/
                   dnorm_mix(data[i], pi, medias, dps),
                 pi * ((data[i] - medias) * 
                         exp(-(data[i] - medias)^2/(2 * dps^2)) / 
                         (sqrt(4 * acos(0)) * dps^3))/
                   dnorm_mix(data[i], pi, medias, dps),
                 pi * (exp(-(data[i] - medias)^2/(2 * dps^2)) * 
                         ((data[i] - medias)^2/(dps^2) - 1) / 
                         (2 *sqrt(4*acos(0))*dps^3))/
                   dnorm_mix(data[i], pi, medias, dps))))
      rownames(si) = c(paste0("pi_", as.character(1:(g-1))), 
                       paste0("mu_", as.character(1:(g))), 
                       paste0("sigma2_", as.character(1:(g))))
      si %*% t(si)
    }
    se = sqrt(diag(solve(Reduce('+', sapply(1:n, si, simplify = FALSE)))))
    class = kmeans(data, centers = medias)$cluster
    if(plot.it){
      output = list(class, pi, medias, dps, se, LF_new, 
                    aic, bic, count, p)
      names(output) = c("classification", "pi_hat", "mu_hat", "sigma_hat",
                        "stde","logLik", "AIC", "BIC", 
                        "EM_iterations", "plot")}
    else{
      output = list(class, pi, medias, dps, se, LF_new, 
                    aic, bic, count)
      names(output) = c("classification", "pi_hat", "mu_hat", "sigma_hat",
                        "stde", "logLik", "AIC", "BIC", 
                        "EM_iterations")
    }
    return(output)
  }
}