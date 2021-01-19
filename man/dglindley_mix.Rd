\name{Density Generalized Lindley Mixture}

\alias{dglindley_mix}

\title{ Generalized Lindley Mixture Density}

\description{
Gives the density value of a univariate Generalized Lindley Mixture.
}

\usage{dglindley_mix(x, pi, alpha, beta, gamma, log = FALSE)}

\arguments{
  \item{x}{ value or vector of values, to apply in the mixture probability density function.}
  \item{pi}{ vector of mixture weights.}
  \item{alpha}{ vector of parameters.}
  \item{beta}{ vector of parameters.}
  \item{gamma}{ vector of parameters.}
  \item{log}{ logical; if TRUE, \eqn{probabilities} p are given as \eqn{log(p)}.}
}

\details{
Each component of the mixture has density
\deqn{f_j(x) = (x^(\alpha_j - 1) * (\alpha_j + \gamma_j*x) * exp(-x/\beta_j))/(\beta_j^(\alpha_j * (\beta_j*\gamma_j + 1) * Gamma(\alpha_j + 1)))}
where \eqn{\alpha_j}, \eqn{\beta_j} and \eqn{\gamma_j} are the parameters of the jth component. j = 1, 2, ..., g.
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
f <- function(x){dglindley_mix(x, pi = c(0.1, 0.9), alpha = c(1, 3), 
                               beta = c(5, 2), gamma = c(1, 2))}
## Plot using graphics package
require(graphics)
curve(f, 0, 25)

require(ggplot2)
## Plot using ggplot2 package
ggplot(data.frame(x=c(0, 25)), aes(x)) +
  stat_function(fun=f, size = 1L, colour = "#0c4c8a") +
  theme_minimal()

## Integrate
integrate(f, 0, "Inf")$value
}
