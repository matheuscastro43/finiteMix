\name{Random Gamma Mixture}

\alias{rgamma_mix}

\title{ Random Generator of a Gamma Mixture}

\description{ Generates random samples of a univariate Gamma Mixture.}

\usage{ rgamma_mix(n, pi, alpha, beta, plot.it = TRUE, empirical = FALSE, 
            col.pop = "red3", col.empirical = "navy", ...) }

\arguments{
  \item{n}{ number of observations.}
  \item{pi}{ vector of mixture weights.}
  \item{alpha}{ shape parameters.}
  \item{beta}{ scale parameters.}
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
data = rgamma_mix(n = 1000, pi = c(0.6, 0.4), alpha = c(2, 4), beta = c(1, 2))

## Sample vector.
data$sample
## The classification of each observation.
data$classification
## The histogram of the sample with population density curve.
data$plot

## Not plotting the graphic.
rgamma_mix(n = 1000, pi = c(0.6, 0.4), alpha = c(2, 4), beta = c(1, 2), plot.it = FALSE)

## Adding the empirical curve to the graphic.
rgamma_mix(n = 1000, pi = c(0.6, 0.4), alpha = c(2, 4), beta = c(1, 2), plot.it = TRUE,
          empirical = TRUE)

## Changing the color of the curves.
rgamma_mix(n = 1000, pi = c(0.6, 0.4), alpha = c(2, 4), beta = c(1, 2), plot.it = TRUE,
          empirical = TRUE, col.pop = "blue", col.empirical = "green")

## Using "...".
rgamma_mix(n = 1000, pi = c(0.6, 0.4), alpha = c(2, 4), beta = c(1, 2), breaks = 300)
}
