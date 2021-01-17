\name{Probabilities Normal Mixture}

\alias{pnorm_mix}

\title{ Distribution Function of a Normal Mixture}

\description{
Gives the distribution function value of a univariate Normal mixture.
}

\usage{
pnorm_mix(q, pi, mean, sd, lower.tail = TRUE)
}

\arguments{
  \item{q}{ value or vector to apply on the distribution function.}
  \item{pi}{ vector of mixture weights.}
  \item{mean}{ vector of mixture means.}
  \item{sd}{ vector of mixture standard desviations.}
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
pnorm_mix(12, pi = c(0.6, 0.4), mean = c(10, 18), sd = c(1, 2))

## A vector
pnorm_mix(c(-100, 10.96684, 100), pi = c(0.6, 0.4), mean = c(10, 18), sd = c(1, 2))


## Survival Function

## A value
pnorm_mix(12, pi = c(0.6, 0.4), mean = c(10, 18), sd = c(1, 2), lower.tail = FALSE)

## A vector
pnorm_mix(c(-100, 10.96684, 100), pi = c(0.6, 0.4), mean = c(10, 18), sd = c(1, 2),
          lower.tail = FALSE)

}
