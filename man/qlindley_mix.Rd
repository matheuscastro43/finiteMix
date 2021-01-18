\name{ Quantile Lindley Mixture}

\alias{qlindley_mix}

\title{ Quantile Function of a Lindley Mixture}

\description{ Gives the quantile function value of a univariate Lindley Mixture.}

\usage{qlindley_mix(p, pi, beta, lower.tail = TRUE, log.p = FALSE) }

\arguments{
  \item{p}{ value or vector to apply on the quantile function.}
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
## Quantile Function

## A value
qlindley_mix(0.59, pi = c(0.6, 0.4), beta = c(1, 5))

## A vector
qlindley_mix(c(0.25, 0.50, 0.75), pi = c(0.6, 0.4), beta = c(1, 5), lower.tail = FALSE)
}
