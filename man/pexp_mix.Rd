\name{Probabilities Exponential Mixture}

\alias{pexp_mix}

\title{ Distribution Function of a Exponential Mixture}

\description{
Gives the distribution function value of a univariate Exponential mixture.
}

\usage{pexp_mix(q, pi, rate, lower.tail = TRUE)}

\arguments{
  \item{q}{ value or vector to apply on the distribution function.}
  \item{pi}{ vector of mixture weights.}
  \item{rate}{ vector of rate parameters.}
  \item{lower.tail}{logical; if TRUE (default). probabilites are \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.}
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
pexp_mix(2, pi = c(0.1, 0.9), rate = c(1, 1/30))

## A vector
pexp_mix(c(0, 17.6336, Inf), pi = c(0.1, 0.9), rate = c(1, 1/30))


## Survival Function

## A value
pexp_mix(2, pi = c(0.8, 0.2), rate = c(1, 1/30), lower.tail = FALSE)

## A vector
pexp_mix(c(0, 10, 15), pi = c(0.1, 0.9), rate = c(1, 1/30), lower.tail = FALSE)
}
