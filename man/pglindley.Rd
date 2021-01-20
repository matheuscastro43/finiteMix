\name{Probabilities Generalized Lindley Distribution}

\alias{pglindley}

\title{ Distribution Function of a Generalized Lindley Distribution}

\description{
Gives the distribution function value of a univariate Generalized Lindley Distribution.
}

\usage{pglindley(q, alpha, beta, gamma, lower.tail = TRUE, log.p = FALSE)}

\arguments{
  \item{q}{ value or vector to apply on the distribution function.}
  \item{alpha}{ value of parameter.}
  \item{beta}{ value of parameter.}
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
pglindley(16, alpha = 2, beta = 4, gamma = 3)

## A vector
pglindley(c(8, 10, 12), alpha = 2, beta = 4, gamma = 3)


## Survival Function

## A value
pglindley(8, alpha = 2, beta = 4, gamma = 3, lower.tail = FALSE)

## A vector
pglindley(c(5, 8, 10), alpha = 2, beta = 4, gamma = 3, lower.tail = FALSE)
}