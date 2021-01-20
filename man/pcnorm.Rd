\name{Probabilities Contaminated Normal Distribution}

\alias{pcnorm}

\title{ Distribution Function of a Contaminated Normal Distribution}

\description{
Gives the distribution function value of a univariate Contaminated Normal Distribution.
}

\usage{pcnorm(q, pi, mean, sd, gamma, lower.tail = TRUE, log.p = FALSE)}

\arguments{
  \item{q}{ value or vector to apply on the distribution function.}
  \item{pi}{ vector of mixture weights.}
  \item{mean}{ value of mixture mean.}
  \item{sd}{ value of parameter.}
  \item{gamma}{ value of parameter.}
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
pcnorm(11.3, pi = c(0.1, 0.9), mean = 10, sd = 2, gamma = 3)

## A vector
pcnorm(c(8, 10, 12), pi = c(0.1, 0.9), mean = 10, sd = 2, gamma = 3)


## Survival Function

## A value
pcnorm(8.7, pi = c(0.8, 0.2), mean = 10, sd = 2, gamma = 3, lower.tail = FALSE)

## A vector
pcnorm(c(6.3, 11, 13.7), pi = c(0.1, 0.9), mean = 10, sd = 2, gamma = 3, lower.tail = FALSE)
}
