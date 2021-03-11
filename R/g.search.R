g.search <- function(data, g = NULL, family = NULL, lim.em = 40, plot.it = TRUE, 
                     col.wss = "navy", col.ic = c("red", "orange"), ...){
  if(is.list(data)){data = data$sample}
  if(length(g) == 0 || length(g) == 1){
    if(is.null(g) || g <= 3){g = 1:3}
    else{g = (g-1):(g+1)}
  }
  wss = sapply(g, function(g){sum(kmeans(data, g)$withinss)})
  pe = noquote(paste0(round(100*(wss/wss[1]), 2), "%"))
  names(wss) = names(pe) = paste("g =", g)
  if(!is.null(family)){
    progress <- function (x, max = length(g)){
      percent <- x / max * 100
      cat(sprintf('\r[%-50s] %d%%',
                  paste(rep('=', percent / 2), collapse = ''),
                  floor(percent)))
      if (x == max)
        cat('\n')
    }
    cat("Estimations (", length(g) ,"): \n", sep = "")
    if(family == "Exponential"){
      aic = bic = NULL
      for(i in g){
        progress(max(c(0, i-1)))
        if(i == 1){
          est = eexp(data, plot.it = FALSE)
        }else{
          suppressMessages(invisible(capture.output(
            est <- eexp_mix(data, g = i, lim.em = lim.em, plot.it = FALSE))))
        }
        aic = c(aic, est$AIC)
        bic = c(bic, est$BIC)
      }
    }
    if(family == "Gamma"){
      aic = bic = NULL
      for(i in g){
        progress(max(c(0, i-1)))
        if(i == 1){
          est = egamma(data, plot.it = FALSE)
        }else{
          suppressMessages(invisible(capture.output(
            est <- egamma_mix(data, g = i, lim.em = lim.em, plot.it = FALSE))))
        }
        aic = c(aic, est$AIC)
        bic = c(bic, est$BIC)
      }
    }
    if(family == "Generalized Lindley"){
      aic = bic = NULL
      for(i in g){
        progress(max(c(0, i-1)))
        if(i == 1){
          est = eglindley(data, plot.it = FALSE)
        }else{
          suppressMessages(invisible(capture.output(
            est <- eglindley_mix(data, g = i, lim.em = lim.em, 
                                 plot.it = FALSE))))
        }
        aic = c(aic, est$AIC)
        bic = c(bic, est$BIC)
      }
    }
    if(family == "Lindley"){
      aic = bic = NULL
      for(i in g){
        progress(max(c(0, i-1)))
        if(i == 1){
          est = elindley(data, plot.it = FALSE)
        }else{
          suppressMessages(invisible(capture.output(
            est <- elindley_mix(data, g = i, lim.em = lim.em, 
                                plot.it = FALSE))))
        }
        aic = c(aic, est$AIC)
        bic = c(bic, est$BIC)
      }
    }
    if(family == "Normal"){
      aic = bic = NULL
      for(i in g){
        progress(max(c(0, i-1)))
        if(i == 1){
          est = enorm(data, plot.it = FALSE)
        }else{
          suppressMessages(invisible(capture.output(
            est <- enorm_mix(data, g = i, lim.em = lim.em, plot.it = FALSE))))
        }
        aic = c(aic, est$AIC)
        bic = c(bic, est$BIC)
      }
    }
    if(family == "Poisson"){
      aic = bic = NULL
      for(i in g){
        progress(max(c(0, i-1)))
        if(i == 1){
          est = epois(data, plot.it = FALSE)
        }else{
          suppressMessages(invisible(capture.output(
            est <- epois_mix(data, g = i, lim.em = lim.em, plot.it = FALSE))))
        }
        aic = c(aic, est$AIC)
        bic = c(bic, est$BIC)
      }
    }
    if(family == "Weibull"){
      aic = bic = NULL
      for(i in g){
        progress(max(c(0, i-1)))
        if(i == 1){
          est = eweibull(data, plot.it = FALSE)
        }else{
          suppressMessages(invisible(capture.output(
            est <- eweibull_mix(data, g = i, lim.em = lim.em, 
                                plot.it = FALSE))))
        }
        aic = c(aic, est$AIC)
        bic = c(bic, est$BIC)
      }
    }
    progress(length(g))
    names(aic) = names(bic) = paste("g =", g)
    aicp = (aic - min(aic))
    aicp = aicp/max(aicp)
    bicp = (bic - min(bic))
    bicp = bicp/max(bicp)
  }
  if(plot.it){
    wssp = wss/wss[1]
    plot(x = g, y = wssp, type = "o", col = col.wss, 
         main = paste0("Total Within-Cluster Sum of Squares", 
                       if(!is.null(family)){" and Infomation Criterions"}), 
         lwd = 4, ylab = "", xaxt = "n", ylim = c(0, 1), ...)
    if(!is.null(family)){
      lines(x = g, aicp, col = col.ic[2], lwd = 5, type = "o")
      lines(x = g, bicp, col = col.ic[1], lwd = 3, type = "o")
      legend("topright", legend=(c("WSS", "AIC/BIC")), 
             fill=c(col.wss, col.ic[1]), 
             border = c(col.wss, col.ic[2]), bty="n")
    }
    else{
      legend("topright", legend="WSS", fill= col.wss, border = col.wss, bty="n")
    }
    axis(1, at = g)
    p <- recordPlot()
    if(!is.null(family)){
      output = list(wss, pe, aic, bic, p)
      names(output) = c("sum","percentageSum", "AIC", "BIC", "plot")
    }else{
      output = list(wss, pe, p)
      names(output) = c("sum","percentageSum", "plot")
    }
  }
  else{
    if(!is.null(family)){
      output = list(wss, pe, aic, bic)
      names(output) = c("sum","percentageSum", "AIC", "BIC")
    }else{
      output = list(wss, pe)
      names(output) = c("sum","percentageSum")
    }
  }
  return(output)
}
