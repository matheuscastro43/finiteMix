\name{ Quantile Lindley Distribution}

\alias{qlindley}

\title{ Quantile Function of a Lindley Distribution}

\description{ Gives the quantile function value of a univariate Lindley Distribution.}

\usage{qlindley(p, beta, lower.tail = TRUE, log.p = FALSE)}

\arguments{
  \item{p}{ value or vector to apply on the quantile function.}
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
## Quantile Function

## A value
qlindley(0.59, beta = 5)

## A vector
qlindley(c(0.25, 0.50, 0.75), beta = 5, lower.tail = FALSE)
}
