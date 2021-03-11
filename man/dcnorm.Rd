\name{Density Contaminated Normal Distribution}

\alias{dcnorm}

\title{ Contaminated Normal Distribution Density}

\description{
Gives the density value of a univariate Contaminated Normal Distribution.
}

\usage{dcnorm(x, pi, mean, sd, gamma, log = FALSE)}

\arguments{
  \item{x}{ value or vector of values, to apply in the mixture probability density function.}
  \item{pi}{ vector of mixture weights.}
  \item{mean}{ value of mixture mean.}
  \item{sd}{ value of parameter.}
  \item{gamma}{ value of parameter.}
  \item{log}{ logical; if TRUE, \eqn{probabilities} p are given as \eqn{log(p)}.}
}

\details{
The mixture has density 
\deqn{f(x) = \pi_1 * N(mean, sd/sqrt(gamma)) + \pi_2 * N(mean, sd)}
where \eqn{N(mean, sd)} is the density function of a Normal Distribution.
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

\author{
CASTRO, M. O.; MONTALVO, G. S. A.}


%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{
f <- function(x){dcnorm(x, pi = c(0.7, 0.3), mean = 10, sd = 1, gamma = 3)}
## Plot using graphics package
require(graphics)
curve(f, 7, 13)

require(ggplot2)
## Plot using ggplot2 package
ggplot(data.frame(x=c(7, 13)), aes(x)) +
  stat_function(fun=f, size = 1L, colour = "#0c4c8a") +
  theme_minimal()

## Integrate
integrate(f, 0, Inf)$value
}
