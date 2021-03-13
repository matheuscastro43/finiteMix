eexp_mix <- function(data, g, lim.em = 100, criteria = "dif.psi", 
                     epsilon = 1e-05, plot.it = TRUE, empirical = FALSE,
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

    L <- function(i){dexp_mix(data[i], pi, rate, log = TRUE)}
    LF <- sum(sapply(1:n, L))

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
        cat("Limit of EM Iterations (", lim.em ,"): \n", sep = "")
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
      pi = pi/sum(pi)
      for(j in 1:g){
        rate[j] <- Wj[j]/(sum(data*Wij[,j]))
      }
      psi_new <- matrix(c(pi, rate), 2, byrow = T)
      LF_new <- sum(sapply(1:n, L))
      
      if(criteria == "dif.lh"){
        crit <- LF_new - LF
        if((abs(crit) < epsilon)){cat("\n"); break}
        LF <- LF_new
      }
      else{
        crit = max(abs(psi - psi_new))
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
    p = prod(dim(matrix(0, 2, 2))) - 1
    aic = 2*p - 2*LF_new
    bic = p*log(n) - 2*LF_new
    if(plot.it == TRUE){
      d.breaks = ceiling(nclass.Sturges(data)*2.5)
      modal = max(hist(data, plot = FALSE, if(any(names(list(...)) == "breaks") == FALSE){
        breaks = d.breaks}, ...)$density)
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
    ordem = order(rate, decreasing = TRUE)
    pi = pi[ordem]
    rate = rate[ordem]
    si = function(i){
      si = t(t(c( (dexp(data[i], rate[-g]) - dexp(data[i], rate[g]))/
                    dexp_mix(data[i], pi, rate),
                  
                  pi * (exp(-rate * data[i]) * (1 - data[i]))/
                    dexp_mix(data[i], pi, rate))))
      rownames(si) = c(paste0("pi_", as.character(1:(g-1))), 
                       paste0("lambda_", as.character(1:(g))))
      si %*% t(si)
    }
    se = sqrt(diag(solve(Reduce('+', sapply(1:n, si, simplify = FALSE)))))
    class = kmeans(data, centers = 1/as.numeric(rate))$cluster
    if(plot.it){
      saida = list(class, pi, as.numeric(rate), se, LF_new, aic,
                   bic, count, p)
      names(saida) = c("classification", "pi_hat", "lambda_hat",
                       "stde", "logLik", "AIC", "BIC", "EM_iterations", "plot")
    }else{
      saida = list(class, pi, as.numeric(rate), se, LF_new, aic,
                   bic, count)
      names(saida) = c("classification", "pi_hat", "lambda_hat", 
                       "stde", "logLik", "AIC", "BIC", "EM_iterations")
    }
    return(saida)
  }
}
