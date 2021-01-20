\name{Probabilities Poisson Mixture}

\alias{ppois_mix}

\title{ Distribution Function of a Poisson Mixture}

\description{
Gives the distribution function value of a univariate Poisson mixture.
}

\usage{
ppois_mix(q, pi, lambda, lower.tail = TRUE, log.p = FALSE)
}

\arguments{
  \item{q}{ value or vector to apply on the distribution function.}
  \item{pi}{ vector of mixture weights.}
  \item{lambda}{ vector of (non-negative) parameters.}
  \item{lower.tail}{logical; if TRUE (default). probabilites are \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.}
  \item{log.p}{logical; if FALSE (default), \eqn{probabilities} are expected in \eqn{p}, otherwise \eqn{log[probabilities]} are expected.}
}

\author{CASTRO, M. O.; MONTALVO, G. S. A.}

\examples{
## Distribution Function

## A value
ppois_mix(2, pi = c(0.1, 0.9), lambda = c(1, 30))

## A vector
ppois_mix(c(0, 17.6336, Inf), pi = c(0.1, 0.9), lambda = c(1, 30))


## Survival Function

## A value
ppois_mix(2, pi = c(0.8, 0.2), lambda = c(1, 30), lower.tail = FALSE)

## A vector
ppois_mix(c(0, 10, 15), pi = c(0.1, 0.9), lambda = c(1, 30), lower.tail = FALSE)
}
