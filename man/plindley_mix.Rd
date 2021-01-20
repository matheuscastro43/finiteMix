\name{Probabilities Lindley Mixture}

\alias{plindley_mix}

\title{ Distribution Function of a Lindley Mixture}

\description{
Gives the distribution function value of a univariate Lindley Mixture.
}

\usage{
plindley_mix(q, pi, beta, lower.tail = TRUE, log.p = FALSE)
}

\arguments{
  \item{q}{ value or vector to apply on the distribution function.}
  \item{pi}{ vector of mixture weights.}
  \item{beta}{ vector of parameters.}
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
plindley_mix(4, pi = c(0.6, 0.4), beta = c(1, 3))

## A vector
plindley_mix(c(4, 7, 11), pi = c(0.6, 0.4), beta = c(1, 3))


## Survival Function

## A value
plindley_mix(2, pi = c(0.6, 0.4), beta = c(1, 3), lower.tail = FALSE)

## A vector
plindley_mix(c(1, 2, 6), pi = c(0.6, 0.4), beta = c(1, 3), lower.tail = FALSE)
}