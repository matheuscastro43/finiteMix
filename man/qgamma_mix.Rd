\name{ Quantile Gamma Mixture}

\alias{qgamma_mix}

\title{ Quantile Function of a Gamma Mixture}

\description{ Gives the quantile function value of a univariate Gamma Mixture.}

\usage{qgamma_mix(p, pi, alpha, beta, lower.tail = TRUE, log.p = FALSE)}

\arguments{
  \item{p}{ value or vector to apply on the quantile function.}
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
## Quantile Function

## A value
qgamma_mix(0.59, pi = c(0.6, 0.4), alpha = c(1, 2), beta = c(1, 5))

## A vector
qgamma_mix(c(0.25, 0.50, 0.75), pi = c(0.6, 0.4), alpha = c(1, 2), beta=c(1, 5))
}
