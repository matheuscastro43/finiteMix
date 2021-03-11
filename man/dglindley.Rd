\name{Density Generalized Lindley Distribution}

\alias{dglindley}

\title{ Generalized Lindley Distribution Density}

\description{
Gives the density value of a univariate Generalized Lindley Distribution.
}

\usage{dglindley(x, alpha, beta, gamma, log = FALSE)}

\arguments{
  \item{x}{ value or vector of values, to apply in the mixture probability density function.}
  \item{alpha}{ value of parameter.}
  \item{beta}{ value of parameter.}
  \item{gamma}{ value of parameter.}
  \item{log}{ logical; if TRUE, \eqn{probabilities} p are given as \eqn{log(p)}.}
}

\details{
The distribution has density 
\deqn{f(x) = (x^(\alpha - 1) * (\alpha + \gamma*x) * exp(-x/\beta))/(\beta^(\alpha * (\beta*\gamma + 1) * Gamma(\alpha + 1))).}
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
f <- function(x){dglindley(x, alpha = 2, beta = 4, gamma = 3)}
## Plot using graphics package
require(graphics)
curve(f, 0, 30)

require(ggplot2)
## Plot using ggplot2 package
ggplot(data.frame(x=c(0, 30)), aes(x)) +
  stat_function(fun=f, size = 1L, colour = "#0c4c8a") +
  theme_minimal()

## Integrate
integrate(f, 0, Inf)$value
}
