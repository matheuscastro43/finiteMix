\name{Density Exponential Mixture}

\alias{dexp_mix}

\title{ Exponential Mixture Density}

\description{
Gives the density value of a univariate Exponential mixture.
}

\usage{dexp_mix(x, pi, rate)}

\arguments{
  \item{x}{ value or vector of values, to apply in the mixture probability density function.}
  \item{pi}{ vector of mixture weights.}
  \item{rate}{ vector of rate parameters.}
}

\details{
Each component of the mixture has density
\deqn{f_j(x) = \lambda_j exp(-\lambda_j * x)}
where \eqn{\lambda_j} is the rate paramter of the jth component. j = 1, 2, ..., g.
}

%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}

%\references{
%% ~put references to the literature/web site here ~
%}

\author{CASTRO, M. O.; MONTALVO, G. S. A.}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{
f <- function(x){dexp_mix(x, pi = c(0.1, 0.9), rate = c(1, 1/30))}
## Plot using graphics package
require(graphics)
curve(f, 0, 40)

require(ggplot2)
## Plot using ggplot2 package
ggplot(data.frame(x=c(0, 40)), aes(x)) +
  stat_function(fun=f, size = 1L, colour = "#0c4c8a") +
  theme_minimal()

## Integrate
integrate(f, 0, "Inf")$value
}
