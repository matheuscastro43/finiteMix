\name{Probabilities Beta Rectangular Distribution}

\alias{pbetar}

\title{ Distribution Function of a Beta Rectangular Distribution}

\description{
Gives the distribution function value of a univariate Beta Rectangular Distribution.
}

\usage{pbetar(q, pi, mu, phi, lower.tail = TRUE, log.p = FALSE)}

\arguments{
  \item{q}{ value or vector to apply on the distribution function.}
  \item{pi}{ vector of mixture weights.}
  \item{mu}{ value of parameter.}
  \item{phi}{ value of parameter.}
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
pbetar(0.3, pi = c(0.1, 0.9), mu = 0.7, phi = 3)

## A vector
pbetar(c(0, 0.6, 1), pi = c(0.1, 0.9), mu = 0.7, phi = 3)


## Survival Function

## A value
pbetar(0.7, pi = c(0.8, 0.2), mu = 0.7, phi = 3, lower.tail = FALSE)

## A vector
pbetar(c(0, 0.10, 0.15), pi = c(0.1, 0.9), mu = 0.7, phi = 3, lower.tail = FALSE)
}
