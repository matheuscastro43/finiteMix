\name{Density Lindley Mixture}

\alias{dlindley_mix}

\title{ Lindley Mixture Density}

\description{
Gives the density value of a univariate Lindley Mixture.
}

\usage{dlindley_mix(x, pi, beta, log = FALSE)}

\arguments{
  \item{x}{ value or vector of values, to apply in the mixture probability density function.}
  \item{pi}{ vector of mixture weights.}
  \item{beta}{ vector of parameters.}
  \item{log}{ logical; if TRUE, \eqn{probabilities} p are given as \eqn{log(p)}.}
}

\details{
Each component of the mixture has density
\deqn{f_j(x) = (x + 1)/(\beta_j * (\beta_j + 1)) * exp(-x/\beta_j)}
where \eqn{\beta_j} is the parameter of the jth component. j = 1, 2, ..., g.
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
f <- function(x){dlindley_mix(x, pi = c(0.1, 0.9), beta = c(5, 2))}
## Plot using graphics package
require(graphics)
curve(f, 0, 15)

require(ggplot2)
## Plot using ggplot2 package
ggplot(data.frame(x=c(0, 15)), aes(x)) +
  stat_function(fun=f, size = 1L, colour = "#0c4c8a") +
  theme_minimal()

## Integrate
integrate(f, 0, "Inf")$value
}
