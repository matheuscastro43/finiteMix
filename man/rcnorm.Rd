\name{Random Contaminated Normal Distribution}

\alias{rnorm}

\title{ Random Generator of a Contaminated Normal Distribution}

\description{ Generates random samples of a univariate Contaminated Normal Distribution.}

\usage{ rcnorm(n, pi, mean, sd, gamma, plot.it = TRUE, empirical = FALSE, col.pop = "red3",
        col.empirical = "navy", ...)}

\arguments{
  \item{n}{ number of observations.}
  \item{pi}{ vector of mixture weights.}
  \item{mean}{ value of mixture mean.}
  \item{sd}{ value of standard deviation parameter.}
  \item{gamma}{ value of gamma parameter.}
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
data = rcnorm(n = 1000, pi = c(0.6, 0.4), mean = 10, sd = 1, gamma = 3)

## Sample vector.
data$sample
## The histogram of the sample with population density curve.
data$plot

## Not plotting the graphic.
rcnorm(n = 1000, pi = c(0.6, 0.4), mean = 10, sd = 1, gamma = 3, plot.it = FALSE)

## Adding the empirical curve to the graphic.
rcnorm(n = 1000, pi = c(0.6, 0.4), mean = 10, sd = 1, gamma = 3, plot.it = TRUE,
          empirical = TRUE)

## Changing the color of the curves.
rcnorm(n = 1000, pi = c(0.6, 0.4), mean = 10, sd = 1, gamma = 3, plot.it = TRUE,
          empirical = TRUE, col.pop = "blue", col.empirical = "green")

## Using "...".
rcnorm(n = 1000, pi = c(0.6, 0.4), mean = 10, sd = 1, gamma = 3, breaks = 300)
}
