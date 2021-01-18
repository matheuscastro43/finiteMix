\name{ Quantile Normal Mixture}

\alias{qnorm_mix}

\title{ Quantile Function of a Normal Mixture}

\description{ Gives the quantile function value of a univariate Normal Mixture.}

\usage{qnorm_mix(p, pi, mean, sd, lower.tail = TRUE, log.p = FALSE)}

\arguments{
  \item{p}{ value or vector to apply on the quantile function.}
  \item{pi}{ vector of mixture weights.}
  \item{mean}{ vector of mixture means.}
  \item{sd}{ vector of mixture standard deviations.}
  \item{lower.tail}{logical; if TRUE (default), \eqn{probabilites} are \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.}
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
qnorm_mix(0.59, pi = c(0.6, 0.4), mean = c(10, 18), sd = c(1, 2))

## A vector
qnorm_mix(c(0.25, 0.50, 0.75), pi = c(0.6, 0.4), mean = c(10, 18), sd = c(1, 2), lower.tail = FALSE)
}
