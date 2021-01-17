\name{Random Normal Mixture}

\alias{rnorm_mix}

\title{ Random Generator of a Normal Mixture}

\description{ Generates random samples of a univariate Normal mixture.}

\usage{ rnorm_mix(n, pi, mean, sd, plot.it = TRUE, empirical = FALSE, col.pop = "red3",
           col.empirical = "navy", ...)}

\arguments{
  \item{n}{ number of observations.}
  \item{pi}{ vector of mixture weights.}
  \item{mean}{ vector of mixture means.}
  \item{sd}{ vector of mixture standard desviations.}
  \item{plot.it}{ logical, TRUE to plot the histogram with estimated distribution curve.}
  \item{empirical}{ logical, TRUE to add the empirical curve ("Kernel Density Estimation") in the plot.}
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
data = rnorm_mix(n = 1000, pi = c(0.6, 0.4), mean = c(10, 18), sd = c(1, 2))

## Sample vector.
data$sample
## The classification of each observation.
data$classification
## The histogram of the sample with population density curve.
data$plot

## Not plotting the graphic.
rnorm_mix(n = 1000, pi = c(0.6, 0.4), mean = c(10, 18), sd = c(1, 2), plot.it = FALSE)

## Adding the empirical curve to the graphic.
rnorm_mix(n = 1000, pi = c(0.6, 0.4), mean = c(10, 18), sd = c(1, 2), plot.it = TRUE,
          empirical = TRUE)

## Changing the color of the curves.
rnorm_mix(n = 1000, pi = c(0.6, 0.4), mean = c(10, 18), sd = c(1, 2), plot.it = TRUE,
          empirical = TRUE, col.pop = "blue", col.empirical = "green")

## Using "...".
rnorm_mix(n = 1000, pi = c(0.6, 0.4), mean = c(10, 18), sd = c(1, 2), breaks = 300)
}
