eexp_mix <- function(data, g, lim.em = 100, criteria = "dif.psi", 
                     plot.it = TRUE, empirical = FALSE,
                     col.estimated = "orange", col.empirical = "navy", ...){
  if((is.numeric(data) || is.numeric(data$sample)) && g == floor(g) && g > 1 &&
     is.logical(plot.it) && is.logical(empirical) &&
     (criteria == "dif.lh" || criteria == "dif.psi")){

    if(is.list(data)){data <- data$sample}
    data <- sort(data)
    n <- length(data)

    k <- kmeans(data, g)
    pi <- table(k$cluster)/n
    rate <- (tapply(data, k$cluster, mean))^(-1)
    psi <- matrix(c(pi, rate), 2, byrow = T)

    L <- function(I) {sum(log(pi[I] * dexp(data[k$cluster == I], rate[I])))}
    LF <- sum(as.numeric(lapply(1:g, L)))

    count = 0
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
          Wij[i,j] <- as.numeric((pi[j]*dexp(data[i], rate[j]))/
                                   sum((pi * dexp(data[i], rate))))
        }
      }
      Wj <- colSums(Wij)
      pi <- 1/n * Wj
      for(j in 1:g){
        rate[j] <- Wj[j]/(sum(data*Wij[,j]))
      }
      psi_new <- matrix(c(pi, rate), 2, byrow = T)
      LF_new <- sum(as.numeric(lapply(1:g, L)))
      
      if(criteria == "dif.lh"){
        crit <- LF_new - LF
        if((abs(crit) < 1*10^(-5)))break;
        LF <- LF_new
      }
      else{
        crit = max(abs(psi - psi_new))
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
      hist(data,freq = F,border = "gray48",
           main = "Sampling distribution of X",xlab = "x",
           ylab = "Density",
           if(any(names(list(...)) == "breaks") == FALSE){breaks = d.breaks}, ...)
      estimada = function(x){dexp_mix(x, pi, rate)}
      curve(estimada, col = col.estimated, lwd = 3, add = T)
      if(empirical){
        lines(density(data),col = col.empirical,lwd = 3)
        legend("topright", legend=(c("Empirical", "Estimated")), fill=c(col.empirical, col.estimated),
               border = c(col.empirical, col.estimated), bty="n")
      }
      else{
        legend("topright", legend = "Estimated", fill = col.estimated, border = col.estimated, bty="n")
      }
      p <- recordPlot()
    }
    ordem = order(rate, decreasing = T)
    class = kmeans(data, centers = 1/as.numeric(rate[ordem]))$cluster
    if(plot.it){
      saida = list(class, pi[ordem], as.numeric(rate[ordem]), count, p)
      names(saida) = c("classification", "pi_hat", "lambda_hat", "EM-interactions", "plot")
    }else{
      saida = list(class, pi[ordem], as.numeric(rate[ordem]), count)
      names(saida) = c("classification", "pi_hat", "lambda_hat", "EM-interactions")
    }
    return(saida)
  }
}
