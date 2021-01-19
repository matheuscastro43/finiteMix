\name{Density Beta Rectangular Distribution}

\alias{dbetar}

\title{ Beta Rectangular Distribution Density}

\description{
Gives the density value of a univariate Beta Rectangular Distribution.
}

\usage{dbetar(x, pi, mu, phi, log = FALSE)}

\arguments{
  \item{x}{ value or vector of values, to apply in the mixture probability density function.}
  \item{pi}{ vector of mixture weights.}
  \item{mu}{ value of parameter.}
  \item{phi}{ value of parameter.}
  \item{log}{ logical; if TRUE, \eqn{probabilities} p are given as \eqn{log(p)}.}
}

\details{
The mixture has density 
\deqn{f(x) = \pi_1 + \pi_2 * (x^(mu*phi - 1) * (1 - x)^((1-mu)*phi)).}
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
f <- function(x){dbetar(x, pi = c(0.5, 0.5), mu = 0.4, phi = 3)}
## Plot using graphics package
require(graphics)
curve(f, 0, 1)

require(ggplot2)
## Plot using ggplot2 package
ggplot(data.frame(x=c(0, 1)), aes(x)) +
  stat_function(fun=f, size = 1L, colour = "#0c4c8a") +
  theme_minimal()

## Integrate
integrate(f, 0, 1)$value
}
