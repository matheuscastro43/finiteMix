\name{Probabilities Generalized Lindley Mixture}

\alias{pglindley_mix}

\title{ Distribution Function of a Generalized Lindley Mixture}

\description{
Gives the distribution function value of a univariate Generalized Lindley Mixture.
}

\usage{
pglindley_mix(q, pi, alpha, beta, gamma, lower.tail = TRUE, log.p = FALSE)
}

\arguments{
  \item{q}{ value or vector to apply on the distribution function.}
  \item{pi}{ vector of mixture weights.}
  \item{alpha}{ vector of parameters.}
  \item{beta}{ vector of parameters.}
  \item{gamma}{ vector of parameters.}
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
pglindley_mix(12, pi = c(0.6, 0.4), alpha = c(2, 4), beta = c(1, 3), 
              gamma = c(2, 3))

## A vector
pglindley_mix(c(-100, 10.96684, 100), pi = c(0.6, 0.4), alpha = c(2, 4), 
              beta = c(1, 3), gamma = c(2, 3))


## Survival Function

## A value
pglindley_mix(12, pi = c(0.6, 0.4), alpha = c(2, 4), beta = c(1, 3), 
              gamma = c(2, 3), lower.tail = FALSE)

## A vector
pglindley_mix(c(5, 10, 16), pi = c(0.6, 0.4), alpha = c(2, 4), beta = c(1, 3), 
              gamma = c(2, 3), lower.tail = FALSE)
}