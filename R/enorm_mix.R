enorm_mix = function(data, g, lim.em = 100, criteria = "dif.psi",
                     plot.it = TRUE, empirical = FALSE, 
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
    L <- function(I){ sum ( log(pi[I] * dnorm(data[which(k$cluster == I)], mean = medias[I], sd = dps[I]) ) ) }
    LF <- sum(as.numeric(lapply(1:g, L))); count = 0
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
      Wij <- matrix(0, nrow = n, ncol = g)
      for(i in 1:n){
        for(j in 1:g){
          Wij[i,j] <- as.numeric((pi[j]*dnorm(data[i], mean = medias[j], sd = dps[j]))/
                                   sum((pi * dnorm(data[i], mean = medias, sd = dps))))
        }
      }
      Wj <- colSums(Wij)
      pi <- 1/n * Wj
      medias <- as.numeric(lapply(1:g, function(j){aux = 0; for(i in 1:n){aux = aux + (data[i]*Wij[i,j])/(Wj[j])};
      return(aux)}))
      dps <- sqrt(as.numeric(lapply(1:g, function(j){aux = 0; for(i in 1:n){aux = aux + ((data[i]-medias[j])^2*Wij[i,j])/(Wj[j])};
      return(aux)})))
      psi_new <- matrix(c(pi, medias, dps), 3, byrow = T)
      LF_new <- sum(as.numeric(lapply(1:g, L)))
      if(criteria == "dif.lh"){
        crit <- LF_new - LF
        if((abs(crit) < 1*10^(-5)))break;
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
        if(crit < 1*10^(-5))break;
        psi <- psi_new
      }
      count = count + 1
      if(count >= lim.em){
        progress(count)
        message("\nLimit of Iterations reached!")
        break
      }
    }
    if(plot.it == TRUE){
      d.breaks <- ceiling(nclass.Sturges(data)*2.5)
      modal <- max(dnorm_mix(medias, pi, medias, dps))
      hist(data,freq = F,border = "gray48",
           main = "Sampling distribution of X", xlab = "x",
           ylab = "Density",
           ylim = c(0, modal),
           if(any(names(list(...)) == "breaks") == FALSE){breaks = d.breaks}, ...)
      
      estimada = function(x){dnorm_mix(x, pi, medias, dps)}
      curve(estimada, col = col.estimated, lwd = 3, add = T)
      if(empirical){
        lines(density(data),col = col.empirical, lwd = 3)
        legend("topright", legend=(c("Empirical", "Estimated")),
               fill=c(col.empirical, col.estimated), border = c(col.empirical, col.estimated),
               bty="n")
      }
      else{
        legend("topright", legend = "Estimated", fill = col.estimated, border = col.estimated,
               bty="n")
      }
      p <- recordPlot()
    }
    ordem = order(medias)
    class = kmeans(data, centers = medias[ordem])$cluster
    if(plot.it){
      output = list(class, pi[ordem], medias[ordem], dps[ordem], count, p)
      names(output) = c("classification", "pi_hat", "mu_hat", "sigma_hat", "EM-interactions", "plot")}
    else{
      output = list(class, pi[ordem], medias[ordem], dps[ordem], count)
      names(output) = c("classification", "pi_hat", "mu_hat", "sigma_hat", "EM-interactions")
    }
    return(output)
  }
}