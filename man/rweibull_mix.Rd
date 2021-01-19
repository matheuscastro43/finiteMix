\name{Random Weibull Mixture}

\alias{rweibull_mix}

\title{ Random Generator of a Weibull Mixture}

\description{ Generates random samples of a univariate Weibull Mixture.}

\usage{ rweibull_mix(n, pi, shape, scale, plot.it = TRUE, empirical = FALSE, 
              col.pop = "red3", col.empirical = "navy", ...) }

\arguments{
  \item{n}{ number of observations.}
  \item{pi}{ vector of mixture weights.}
  \item{shape}{ shape parameters.}
  \item{scale}{ scale parameters.}
  \item{plot.it}{ logical; if TRUE (default), plot the histogram with estimated distribution curve.}
  \item{empirical}{ logical; if TRUE, add the empirical curve ("Kernel Density Estimation") in the plot.}
  \item{col.pop}{ a colour to be used in the curve of population density.}
  \item{col.empirical}{ a colour to be used in the curve of empirical density.}
  \item{...}{ further arguments and graphical parameters passed to \link[graphics]{hist}.}
}

%\details{}
%\value{}
%\references{}

\author{CASTRO, M. O.; MONTALVO, G. S. A.}

%\note{}
%\seealso{}

\examples{
## Generate a sample.
data = rweibull_mix(n = 1000, pi = c(0.6, 0.4), shape = c(2, 4), scale = c(1, 2))

## Sample vector.
data$sample
## The classification of each observation.
data$classification
## The histogram of the sample with population density curve.
data$plot

## Not plotting the graphic.
rweibull_mix(n = 1000, pi = c(0.6, 0.4), shape = c(2, 4), scale = c(1, 2), plot.it = FALSE)

## Adding the empirical curve to the graphic.
rweibull_mix(n = 1000, pi = c(0.6, 0.4), shape = c(2, 4), scale = c(1, 2), plot.it = TRUE,
          empirical = TRUE)

## Changing the color of the curves.
rweibull_mix(n = 1000, pi = c(0.6, 0.4), shape = c(2, 4), scale = c(1, 2), plot.it = TRUE,
          empirical = TRUE, col.pop = "blue", col.empirical = "green")

## Using "...".
rweibull_mix(n = 1000, pi = c(0.6, 0.4), shape = c(2, 4), scale = c(1, 2), breaks = 300)
}
