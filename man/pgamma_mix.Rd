\name{Probabilities Gamma Mixture}

\alias{pgamma_mix}

\title{ Distribution Function of a Gamma Mixture}

\description{
Gives the distribution function value of a univariate Gamma Mixture.
}

\usage{
pgamma_mix(q, pi, alpha, beta, lower.tail = TRUE, log.p = FALSE)
}

\arguments{
  \item{q}{ value or vector to apply on the distribution function.}
  \item{pi}{ vector of mixture weights.}
  \item{alpha}{ shape parameters.}
  \item{beta}{ scale parameters.}
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
## Distribution Function

## A value
pgamma_mix(12, pi = c(0.6, 0.4), alpha = c(2, 4), beta = c(1, 3))

## A vector
pgamma_mix(c(-100, 10.96684, 100), pi = c(0.6, 0.4), alpha = c(2, 4), beta = c(1, 3))


## Survival Function

## A value
pgamma_mix(12, pi = c(0.6, 0.4), alpha = c(2, 4), beta = c(1, 3), lower.tail = FALSE)

## A vector
pgamma_mix(c(5, 10, 16), pi = c(0.6, 0.4), alpha = c(2, 4), beta = c(1, 3),
          lower.tail = FALSE)
}