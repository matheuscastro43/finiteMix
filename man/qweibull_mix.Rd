\name{ Quantile Weibull Mixture}

\alias{qweibull_mix}

\title{ Quantile Function of a Weibull Mixture}

\description{ Gives the quantile function value of a univariate Weibull Mixture.}

\usage{qweibull_mix(p, pi, shape, scale, lower.tail = TRUE, log.p = FALSE)}

\arguments{
  \item{p}{ value or vector to apply on the quantile function.}
  \item{pi}{ vector of mixture weights.}
  \item{shape}{ vector of shape parameters.}
  \item{scale}{ vector of scale parameters.}
  \item{lower.tail}{logical; if TRUE (default), \eqn{probabilites} are \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.}
  \item{log.p}{logical; if FALSE (default), \eqn{probabilities} are expected in \eqn{p}, otherwise \eqn{log[probabilities]} are expected.}
}

%\details{}
%\value{}
%\references{}
\author{CASTRO, M. O.; MONTALVO, G. S. A.}
%\note{}
%\seealso{}
\examples{
## Quantile Function

## A value
qweibull_mix(0.59, pi = c(0.6, 0.4), shape = c(10, 18), scale = c(1, 2))

## A vector
qweibull_mix(c(0.25, 0.50, 0.75), pi = c(0.6, 0.4), shape = c(10, 18), scale = c(1, 2))
}
