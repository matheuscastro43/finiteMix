\name{Mode Generalized Lindley Mixture}

\alias{moglindley_mix}

\title{ Gemeralized Lindley Mixture Mode}

\description{
Gives the mode value of a univariate Generalized Lindley Mixture.
}

\usage{moglindley_mix(pi, alpha, beta, gamma)}

\arguments{
  \item{pi}{ vector of mixture weights.}
  \item{alpha}{ vector of parameters.}
  \item{beta}{ vector of parameters.}
  \item{gamma}{ vector of parameters.}
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
f <- function(x){dglindley_mix(x, c(.5, .5), c(1, 2), c(2, 2), c(1,3))}
## Plot using graphics package
require(graphics)
curve(f, 0, 6)
mode <- moglindley_mix(c(.5, .5), c(1, 2), c(2, 2), c(1,3))
lines(c(mode, mode), c(0, f(mode)), lwd = 2)
}
