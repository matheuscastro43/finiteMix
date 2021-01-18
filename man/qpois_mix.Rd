\name{ Quantile Poisson Mixture}

\alias{qpois_mix}

\title{ Quantile Function of a Poisson Mixture}

\description{ Gives the quantile function value of a univariate Poisson Mixture.}

\usage{qpois_mix(p, pi, lambda, lower.tail = TRUE, log.p = FALSE)}

\arguments{
  \item{p}{ value or vector to apply on the quantile function.}
  \item{pi}{ vector of mixture weights.}
  \item{lambda}{ vector of (non-negative) parameters.}
  \item{lower.tail}{logical; if TRUE (default). probabilites are \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.}
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
qpois_mix(0.59, pi = c(0.6, 0.4), lambda = c(10, 18))
qpois_mix(0.59, pi = c(0.6, 0.4), lambda = c(10, 18), lower.tail = FALSE)

## A vector
qpois_mix(c(0.25, 0.50, 0.75), pi = c(0.6, 0.4), lambda = c(10, 18))
qpois_mix(c(0.25, 0.50, 0.75), pi = c(0.6, 0.4), lambda = c(10, 18), lower.tail = FALSE)
}
