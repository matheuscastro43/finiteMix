\name{Probabilities Lindley Distribution}

\alias{plindley}

\title{ Distribution Function of a Generalized Distribution}

\description{
Gives the distribution function value of a univariate Lindley Distribution.
}

\usage{plindley(q, beta, lower.tail = TRUE, log.p = FALSE)}

\arguments{
  \item{q}{ value or vector to apply on the distribution function.}
  \item{beta}{ value of parameter.}
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
plindley(16, beta = 4)

## A vector
plindley(c(8, 10, 12), beta = 4)


## Survival Function

## A value
plindley(8, beta = 4, lower.tail = FALSE)

## A vector
plindley(c(5, 8, 10), beta = 4, lower.tail = FALSE)
}