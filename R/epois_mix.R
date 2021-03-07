epois_mix = function(data, g, lim.em = 100, criteria = "dif.psi", 
                     plot.it = TRUE, empirical = FALSE, 
                     col.estimated = "orange", col.empirical = "navy", ...){
  if((is.numeric(data) || is.numeric(data$sample)) && is.logical(plot.it) &&
     is.logical(empirical) && (criteria == "dif.lh" || criteria == "dif.psi")){
    
    if(is.list(data)){data <- data$sample}
    data <- sort(data)
    n <- length(data)
    
    k <- kmeans(-data, g)
    pi <- table(k$cluster)/n
    medias <- tapply(data, k$cluster, mean)
    variancias <- tapply(data, k$cluster, var)
    lambda <- apply(rbind(medias, variancias), 2, mean)
    psi <- matrix(c(pi, lambda), 2, byrow = T)
    
    L <- function(I){ sum ( log(pi[I] * dpois(data[which(k$cluster == I)],
                                              lambda = lambda[I]) ) ) }
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
          Wij[i,j] <- as.numeric((pi[j]*dpois(data[i], lambda = lambda[j]))/
                                   sum((pi * dpois(data[i], lambda = lambda))))
        }
      }
      Wj <- colSums(Wij)
      lambda_j <- function(j){aux = 0; for(i in 1:n){aux = aux + (data[i]*Wij[i,j])/(Wj[j])};
      return(aux)}
      pi <- 1/n * Wj
      lambda <- as.numeric(lapply(1:g, lambda_j))
      psi_new <- matrix(c(pi, lambda), 2, byrow = T)
      LF_new <- sum(as.numeric(lapply(1:g, L)))
      if(criteria == "dif.lh"){
        crit <- LF_new - LF
        if((abs(crit) < 1*10^(-5))){cat("\n"); break}
        LF <- LF_new
      }
      else{
        crit = max(abs(psi - psi_new))
        if(any(is.na(crit))){
          k <- kmeans(-data, g)
          pi <- table(k$cluster)/n
          medias <- tapply(data, k$cluster, mean)
          variancias <- tapply(data, k$cluster, var)
          lambda <- apply(rbind(medias, variancias), 2, mean)
          psi <- matrix(c(pi, lambda), 2, byrow = T)
          next
        }
        if(crit < 1*10^(-5)){cat("\n"); break}
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
      modal = max(max(dpois_mix(floor(lambda), pi, lambda)),
                  hist(data, plot = FALSE, if(any(names(list(...)) == "breaks") == FALSE){
                    breaks = d.breaks}, ...)$density)
      hist(data,freq = F,border = "gray48",
           main = "Sampling distribution of X",xlab = "x",
           ylab = "Density",
           if(any(names(list(...)) == "breaks") == FALSE){breaks = d.breaks}, ...)
      estimada = function(x){dpois_mix(x, pi, lambda)}
      x0 <- min(data)
      x1 <- max(data)
      par(new = T)
      plot(x0:x1, lapply(x0:x1, estimada), type = "h", col = col.estimated,
           lwd = 3, main = "",
           xlab = "", ylab = "", axes = F)
      if(empirical){
        lines(density(data),col = col.empirical,lwd = 3)
        legend("topright", legend=(c("Empirical", "Estimated")),
               fill=c(col.empirical, col.estimated),
               border = c(col.empirical, col.estimated), bty="n")
      }
      else{
        legend("topright", legend = "Estimated", fill = col.estimated,
               border = col.estimated, bty="n")
      }
      p <- recordPlot()
    }
    ordem = order(medias)
    class = kmeans(data, centers = lambda[ordem])$cluster
    if(plot.it){
      saida = list(class, pi[ordem], lambda[ordem], LF_new, count, p)
      names(saida) = c("classification" ,"pi_hat", "lambda_hat", 
                       "logLik", "EM-interactions", "plot")
    }else{
      saida = list(class, pi[ordem], lambda[ordem], LF_new, count)
      names(saida) = c("classification" ,"pi_hat", "lambda_hat",
                       "logLik", "EM-interactions")
    }
    return(saida)
  }
}
