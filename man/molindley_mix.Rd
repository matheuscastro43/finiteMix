\name{Mode Lindley Mixture}

\alias{molindley_mix}

\title{ Lindley Mixture Mode}

\description{
Gives the mode value of a univariate Lindley Mixture.
}

\usage{molindley_mix(pi, beta)}

\arguments{
  \item{pi}{ vector of mixture weights.}
  \item{beta}{ vector of parameters.}
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
f <- function(x){dlindley_mix(x, c(.5, .5), c(2, 6))}
## Plot using graphics package
require(graphics)
curve(f, 0, 6)
mode <- molindley_mix(c(.5, .5), c(2, 6))
lines(c(mode, mode), c(0, f(mode)), lwd = 2)
}
