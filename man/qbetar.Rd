\name{ Quantile Beta Rectangular Distribution}

\alias{qbetar}

\title{ Quantile Function of a Beta Rectangular Distribution}

\description{ Gives the quantile function value of a univariate Beta Rectangular Distribution.}

\usage{qbetar(p, pi, mu, phi, lower.tail = TRUE, log.p = FALSE)}

\arguments{
  \item{p}{ value or vector to apply on the quantile function.}
  \item{pi}{ vector of mixture weights.}
  \item{mu}{ value of parameter.}
  \item{phi}{ value of parameter.}
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
qbetar(0.59, pi = c(0.6, 0.4), mu = .4, phi = 3)

## A vector
qbetar(c(0.25, 0.50, 0.75), pi = c(0.6, 0.4), mu = .4, phi = 3, lower.tail = FALSE)
}
