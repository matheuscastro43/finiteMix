\name{ Quantile Exponential Mixture}

\alias{qexp_mix}

\title{ Quantile Function of a Exponential Mixture}

\description{ Gives the quantile function value of a univariate Exponential mixture.}

\usage{qexp_mix(p, pi, rate, lower.tail = TRUE, log.p = FALSE)}

\arguments{
  \item{p}{ value or vector to apply on the quantile function.}
  \item{pi}{ vector of mixture weights.}
  \item{rate}{ vector of rate parameters.}
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
qexp_mix(0.59, pi = c(0.6, 0.4), rate = c(1, 5))

## A vector
qexp_mix(c(0.25, 0.50, 0.75), pi = c(0.6, 0.4), rate = c(1, 5))
}
