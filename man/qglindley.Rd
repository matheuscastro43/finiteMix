\name{ Quantile Generalized Lindley Distribution}

\alias{qglindley}

\title{ Quantile Function of a Generalized Lindley Distribution}

\description{ Gives the quantile function value of a univariate Generalized Lindley Distribution.}

\usage{qglindley(p, alpha, beta, gamma, lower.tail = TRUE, log.p = FALSE)}

\arguments{
  \item{p}{ value or vector to apply on the quantile function.}
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
## Quantile Function

## A value
qglindley(0.59, alpha = 2, beta = 5, gamma = 3)

## A vector
qglindley(c(0.50, 0.75), alpha = 2, beta = 5, gamma = 3)
}
