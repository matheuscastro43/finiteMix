\name{Random Beta Rectangular Distribution}

\alias{rbetar}

\title{ Random Generator of a Beta Rectangular Distribution}

\description{ Generates random samples of a univariate Beta Rectangular Distribution.}

\usage{ rbetar(n, pi, mu, phi, plot.it = TRUE, empirical = FALSE, col.pop = "red3",
        col.empirical = "navy", ...)}

\arguments{
  \item{n}{ number of observations.}
  \item{pi}{ vector of mixture weights.}
  \item{mu}{ value of parameter.}
  \item{phi}{ value of parameter.}
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
data = rbetar(n = 1000, pi = c(0.6, 0.4), mu = 0.4, phi = 3)

## Sample vector.
data$sample
## The histogram of the sample with population density curve.
data$plot

## Not plotting the graphic.
rbetar(n = 1000, pi = c(0.6, 0.4), mu = 0.4, phi = 3, plot.it = FALSE)

## Adding the empirical curve to the graphic.
rbetar(n = 1000, pi = c(0.6, 0.4), mu = 0.4, phi = 3, plot.it = TRUE,
       empirical = TRUE)

## Changing the color of the curves.
rbetar(n = 1000, pi = c(0.6, 0.4), mu = 0.4, phi = 3, plot.it = TRUE,
       empirical = TRUE, col.pop = "blue", col.empirical = "green")

## Using "...".
rbetar(n = 1000, pi = c(0.6, 0.4), mu = 0.4, phi = 3, breaks = 300)
}
