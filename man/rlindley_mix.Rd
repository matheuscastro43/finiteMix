\name{Random Lindley Mixture}

\alias{rlindley_mix}

\title{ Random Generator of a Lindley Mixture}

\description{ Generates random samples of a univariate Lindley Mixture.}

\usage{ rlindley_mix(n, pi, beta, plot.it = TRUE, empirical = FALSE, col.pop = "red3",
              col.empirical = "navy", ...)}

\arguments{
  \item{n}{ number of observations.}
  \item{pi}{ vector of mixture weights.}
  \item{beta}{ vector of parameters.}
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
data = rlindley_mix(n = 1000, pi = c(0.6, 0.4), beta = c(2, 3))

## Sample vector.
data$sample
## The histogram of the sample with population density curve.
data$plot

## Not plotting the graphic.
rlindley_mix(n = 1000, pi = c(0.6, 0.4), beta = c(2, 3), plot.it = FALSE)

## Adding the empirical curve to the graphic.
rlindley_mix(n = 1000, pi = c(0.6, 0.4), beta = c(2, 3), plot.it = TRUE, 
             empirical = TRUE)

## Changing the color of the curves.
rlindley_mix(n = 1000, pi = c(0.6, 0.4), beta = c(2, 3), plot.it = TRUE,
             empirical = TRUE, col.pop = "blue",col.empirical = "green")

## Using "...".
rlindley_mix(n = 1000, pi = c(0.6, 0.4), beta = c(2, 3), breaks = 300)
}