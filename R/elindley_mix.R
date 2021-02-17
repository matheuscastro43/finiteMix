elindley_mix = function(data, g, lim.em = 100, criteria = "dif.psi", plot.it = 
                          TRUE, empirical = FALSE, col.estimated = "orange", 
                        col.empirical = "navy", ...){
  if((is.numeric(data) || is.numeric(data$sample)) && g == floor(g) && g > 1 &&
     is.logical(plot.it) && is.logical(empirical) &&
     (criteria == "dif.lh" || criteria == "dif.psi")){
    if(is.list(data)){data = data$sample}
    data = sort(data)
    n = length(data)
    
    psi = matrix(0, ncol = g, nrow = 2)
    k = kmeans(data, g)
    for(j in 1:g){
      est = elindley(data[k$cluster == j], plot.it = F)
      psi[2, j] = c(est$beta_hat)
    }
    betas = psi[2,]
    medias = (betas * (1 + 2*betas))/(1 + betas)
    k = kmeans(data, g, centers = medias)
    psi[1, ] = pi = table(k$cluster)/n
    for(j in 1:g){
      est = elindley(data[k$cluster == j], plot.it = F)
      psi[2, j] = c(est$beta_hat)
    }
    betas = psi[2,]
    
    L = function(I){sum(log(pi[I] * dlindley(data[which(k$cluster == I)], 
                                             beta = betas[I])))}
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
          Wij[i,j] = as.numeric((pi[j]*dlindley(data[i], beta = betas[j]))/
                                  sum((pi * dlindley(data[i], beta = betas))))
        }
      }
      Wj <- colSums(Wij)
      pi <- 1/n * Wj
      
      Q = function(param){
        betast = param
        
        Qj = function(j){
          aux = 0
          for(i in 1:n){
            aux = aux + Wij[i,j] * log(pi[j] * dlindley(data[i], betast[j]))
          }
          return(aux)
        }
        q = sum(sapply(1:g, Qj))
        if(q == -Inf) return(.Machine$double.xmax/1e+08)
        if(q == Inf) return(-.Machine$double.xmax/1e+08)
        return(-q)
      }
      grQ = function(param){
        betast = param
        
        gr = rep(0, g)
        for(i in 1:n){
          gr = gr + Wij[i, ] * data[i]/(betast^2) - Wij[i, ]/betast - 
            Wij[i, j]/(betast + 1)
        }
        return(-gr)
      }
      
      estim = optim(par = betas, fn = Q, method = "L-BFGS-B", lower = 1e-04,
                    upper = Inf, gr = grQ)
      betas = estim$par
      
      psi_new = matrix(c(pi, betas), 2, byrow = T)
      LF_new = sum(sapply(1:g, L))
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
            est = elindley(data[k$cluster == j], plot.it = F)
            psi[2, j] = c(est$beta_hat)
          }
          betas = psi[2,]
          medias = (betas * (1 + 2*betas))/(1 + betas)
          k = kmeans(data, g, centers = medias)
          psi[1, ] = pi = table(k$cluster)/n
          for(j in 1:g){
            est = elindley(data[k$cluster == j], plot.it = F)
            psi[2, j] = c(est$beta_hat)
          }
          betas = psi[2,]
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
    modal = 0
    for(i in 1:g){
      modal[i] = max(dlindley_mix(c(0, betas[i]), pi, betas))
      if(modal[i] > 10* dlindley_mix(1, pi, betas)){
        modal[i] = dlindley_mix(1, pi, betas)
      }
    }
    modal = min(c(1, max(modal, hist(data, if(any(names(list(...)) == 
                                                  "breaks") == FALSE){
      breaks = d.breaks}, ...)$density)))
    hist(data, freq = F, border = "gray48",
         main = "Sampling distribution of X", xlab = "x",
         ylab = "Density",
         ylim = c(0, modal),
         if(any(names(list(...)) == "breaks") == FALSE){breaks = d.breaks}, ...)
    
    estimada = function(x){dlindley_mix(x, pi, betas)}
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
  medias = (betas * (1 + 2*betas))/(1 + betas)
  ordem = order(medias)
  class = kmeans(data, centers = medias[ordem])$cluster
  if(plot.it){
    output = list(class, pi[ordem], betas[ordem], count, p)
    names(output) = c("classification", "pi_hat", "beta_hat", 
                      "EM-interactions", "plot")}
  else{
    output = list(class, pi[ordem], betas[ordem], count)
    names(output) = c("classification", "pi_hat", "beta_hat", 
                      "EM-interactions")
  }
  return(output)
}
