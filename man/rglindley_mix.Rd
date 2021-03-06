\name{Random Generalized Lindley Mixture}

\alias{rglindley_mix}

\title{ Random Generator of a Generalized Lindley Mixture}

\description{ Generates random samples of a univariate Generalizeed Lindley Mixture.}

\usage{ rglindley_mix(n, pi, alpha, beta, gamma, plot.it = TRUE, empirical = FALSE, 
               col.pop = "red3", col.empirical = "navy", ...)}

\arguments{
  \item{n}{ number of observations.}
  \item{pi}{ vector of mixture weights.}
  \item{alpha}{ vector of parameters.}
  \item{beta}{ vector of parameters.}
  \item{gamma}{ vector of parameters.}
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
data = rglindley_mix(n = 1000, pi = c(0.6, 0.4), alpha = c(2, 6),
                     beta = c(2, 3), gamma = c(4, 2))

## Sample vector.
data$sample
## The histogram of the sample with population density curve.
data$plot

## Not plotting the graphic.
rglindley_mix(n = 1000, pi = c(0.6, 0.4), alpha = c(2, 6), beta = c(2, 3),
              gamma = c(4, 2), plot.it = FALSE)

## Adding the empirical curve to the graphic.
rglindley_mix(n = 1000, pi = c(0.6, 0.4), alpha = c(2, 6), beta = c(2, 3),
              gamma = c(4, 2), plot.it = TRUE,
              empirical = TRUE)

## Changing the color of the curves.
rglindley_mix(n = 1000, pi = c(0.6, 0.4), alpha = c(2, 6), beta = c(2, 3),
              gamma = c(4, 2), plot.it = TRUE,empirical = TRUE, 
              col.pop = "blue", col.empirical = "green")

## Using "...".
rglindley_mix(n = 1000, pi = c(0.6, 0.4), alpha = c(2, 6), beta = c(2, 3), 
              gamma = c(4, 2), breaks = 300)
}