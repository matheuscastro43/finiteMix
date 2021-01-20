\name{Probabilities Weibull Mixture}

\alias{pweibull_mix}

\title{ Distribution Function of a Weibull Mixture}

\description{
Gives the distribution function value of a univariate Weibull Mixture.
}

\usage{
pweibull_mix(q, pi, shape, scale, lower.tail = TRUE, log.p = FALSE)
}

\arguments{
  \item{q}{ value or vector to apply on the distribution function.}
  \item{pi}{ vector of mixture weights.}
  \item{shape}{ shape parameters.}
  \item{scale}{ scale parameters.}
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
pweibull_mix(4, pi = c(0.6, 0.4), shape = c(1, 3), scale = c(2, 4))

## A vector
pweibull_mix(c(4, 7, 11), pi = c(0.6, 0.4), shape = c(1, 3), scale = c(2, 4))


## Survival Function

## A value
pweibull_mix(2, pi = c(0.6, 0.4), shape = c(1, 3), scale = c(2, 4),
             lower.tail = FALSE)

## A vector
pweibull_mix(c(1, 2, 6), pi = c(0.6, 0.4), shape = c(1, 3), scale = c(2, 4),
             lower.tail = FALSE)
}