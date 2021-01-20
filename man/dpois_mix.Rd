\name{Probability Mass Poisson Mixture}

\alias{dpois_mix}

\title{ Poisson Mixture Probability Mass}

\description{
Gives the probability mass value of a univariate Poisson Mixture.
}

\usage{dpois_mix(x, pi, lambda, log = FALSE)}

\arguments{
  \item{x}{ value or vector of values, to apply in the mixture probability mass function.}
  \item{pi}{ vector of mixture weights.}
  \item{lambda}{vector os parameters (means).}
  \item{log}{ logical; if TRUE, \eqn{probabilities} p are given as \eqn{log(p)}.}
}

\details{
Each component of the mixture has density
\deqn{f_j(x) = (exp(- \lambda_j) * \lambda_j^(x))) / (x!)}
where \eqn{\lambda_j} is the mean and the variance of the jth component. j = 1, 2, ..., g.
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

%\author{
%%  ~~who you are~~
%}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{
f <- function(x){dpois_mix(x, c(0.5, 0.5), c(3, 8))}
## Plot using graphics package
require(graphics)
plot.function(f, 0, 10, type = "h")

require(ggplot2)
## Plot using ggplot2 package
x <- seq(0, 20, 1)
aux <- data.frame(x, y = f(x))

ggplot(aux, aes(x=x)) +
  geom_point( aes(y=y) ) +
  geom_linerange(aes(ymax=y), ymin=0) +
  theme_minimal()

## Sum
sum(f(0:10000))
}