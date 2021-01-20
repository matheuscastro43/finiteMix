rexp_mix <- function(n, pi, rate, plot.it = TRUE, empirical = FALSE, col.pop = "red3",
                     col.empirical = "navy", ...){
  g <- length(pi)
  if(n == floor(n) && sum(pi) == 1 && min(c(pi, rate, n)) > 0 && length(rate) == g){

    z <- rmultinom(n, 1, pi)
    aux <- rowSums(z)

    sample <- NULL
    for(j in 1:g){
      sample <- c(sample, rexp(aux[j], rate[j]))
    }
    if(plot.it){
      d.breaks <- ceiling(nclass.Sturges(sample)*2.5)
      hist(sample,freq = F,border = "gray48",
           main = "Sampling distribution of X",xlab = "x",
           ylab = "Density",
           if(any(names(list(...)) == "breaks") == FALSE){breaks = d.breaks}, ...)
      pop = function(x){dexp_mix(x, pi, rate)}
      curve(pop, col = col.pop, lwd = 3, add = T)
      if(empirical){
        lines(density(sample),col = col.empirical,lwd = 3)
        legend("topright", legend=(c("Population", "Empirical")),
               fill=c(col.pop, col.empirical),border = c(col.pop, col.empirical), bty="n")
      }
      else{
        legend("topright", legend=(c("Population")),
               fill=c(col.pop),border = c(col.pop), bty="n")
      }
      p <- recordPlot()
    }
    ord <- order(sample)
    sample <- cbind(sample, rep(1:g, aux))
    sample <- sample[ord,]
    if(plot.it){
      output = list(sample[,1], g, pi, rate, sample[,2], p)
      names(output) = c("sample", "g", "pi", "lambda", "classification", "plot")
    }
    else{
      output = list(sample[,1], g, pi, rate, sample[,2])
      names(output) = c("sample", "g", "pi", "lambda", "classification")
    }
    return(output)}
  else stop("The parametric space must be respected.")
}
