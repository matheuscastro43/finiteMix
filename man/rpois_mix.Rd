\name{Random Poisson Mixture}

\alias{rpois_mix}

\title{ Random Generator of a Poisson Mixture}

\description{ Generates random samples of a univariate Poisson Mixture.}

\usage{ rpois_mix(n, pi, lambda, plot.it = TRUE, empirical = FALSE, col.pop = "red3",
                     col.empirical = "navy", ...)}

\arguments{
  \item{n}{ number of observations.}
  \item{pi}{ vector of mixture weights.}
  \item{lambda}{ vector of (non-negative) parameters.}
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
data = rpois_mix(n = 1000, pi = c(0.6, 0.4), lambda = c(10, 15))

## Sample vector.
data$sample
## The classification of each observation.
data$classification
## The histogram of the sample with population density curve.
data$plot

## Not plotting the graphic.
rpois_mix(n = 1000, pi = c(0.6, 0.4), lambda = c(10, 25), plot.it = FALSE)

## Adding the empirical curve to the graphic.
rpois_mix(n = 1000, pi = c(0.6, 0.4), lambda = c(10, 25), plot.it = TRUE,
          empirical = TRUE)

## Changing the color of the curves.
rpois_mix(n = 1000, pi = c(0.6, 0.4), lambda = c(10, 25), plot.it = TRUE,
          empirical = TRUE, col.pop = "blue", col.empirical = "green")

## Using "...".
rpois_mix(n = 1000, pi = c(0.6, 0.4), lambda = c(10, 25), breaks = 300)
}
