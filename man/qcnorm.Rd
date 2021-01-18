\name{ Quantile Contaminated Normal Distribution}

\alias{qcnorm}

\title{ Quantile Function of a Contaminated Normal Distribution}

\description{ Gives the quantile function value of a univariate Contaminated Normal Distribution.}

\usage{qcnorm(p, pi, mean, sd, gamma, lower.tail = TRUE, log.p = FALSE)}

\arguments{
  \item{p}{ value or vector to apply on the quantile function.}
  \item{pi}{ vector of mixture weights.}
  \item{mean}{ value of mixture mean.}
  \item{sd}{ value of standard deviation parameter.}
  \item{gamma}{ value of gamma parameter.}
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
qcnorm(0.59, pi = c(0.6, 0.4), mean = 10, sd = 2, gamma = 3)

## A vector
qcnorm(c(0.25, 0.50, 0.75), pi = c(0.6, 0.4), mean = 10, sd = 2, gamma = 3, lower.tail = FALSE)
}
