\name{Random Exponential Mixture}

\alias{rexp_mix}

\title{ Random Generator of a Exponential Mixture}

\description{ Generates random samples of a univariate Exponential Mixture.}

\usage{ rexp_mix(n, pi, rate, plot.it = TRUE, empirical = FALSE, col.pop = "red3",
         col.empirical = "navy", ...)}

\arguments{
  \item{n}{ number of observations.}
  \item{pi}{ vector of mixture weights.}
  \item{rate}{ vector of rate parameters.}
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
data = rexp_mix(n = 1000, pi = c(0.1, 0.9), rate = c(1, 1/1000))

## Sample vector.
data$sample
## The classification of each observation.
data$classification
## The histogram of the sample with population density curve.
data$plot

## Not plotting the graphic.
rexp_mix(n = 1000, pi = c(0.1, 0.9), rate = c(1, 1/1000), plot.it = FALSE)

## Adding the empirical curve to the graphic.
rexp_mix(n = 1000, pi = c(0.1, 0.9), rate = c(1, 1/1000), empirical = TRUE)

## Changing the color of the curves.
rexp_mix(n = 1000, pi = c(0.1, 0.9), rate = c(1, 1/1000), empirical = TRUE,
         col.pop = "blue", col.empirical = "green")

## Using "...".
rexp_mix(n = 1000, pi = c(0.1, 0.9), rate = c(1, 1/1000), breaks = 300)
}
