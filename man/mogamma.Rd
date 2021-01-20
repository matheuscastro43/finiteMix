\name{Mode Gamma Distribution}

\alias{mogamma}

\title{ Gamma Distribution Mode}

\description{
Gives the mode value of a univariate Gamma Distribution.
}

\usage{mogamma(alpha, beta)}

\arguments{
  \item{alpha}{ shape parameter.}
  \item{beta}{ scale parameter.}
}

\details{
The distribution has density
\deqn{f(x) = 1/(\beta^(\alpha) Gamma(\alpha)) x^(\alpha-1) e^-(x/\beta)}
where \eqn{\alpha} and \eqn{\beta}are the parameters.
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
f <- function(x){dgamma(x, 2, 2)}
## Plot using graphics package
require(graphics)
curve(f, 0, 2)
mode <- mogamma(2, 2)
lines(c(mode, mode), c(0, f(mode)), lwd = 2)
}
