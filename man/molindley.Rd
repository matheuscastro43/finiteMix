\name{Mode Lindley Distribution}

\alias{molindley}

\title{ Lindley Distribution Mode}

\description{
Gives the mode value of a univariate Lindley Distribution.
}

\usage{molindley(beta)}

\arguments{
  \item{beta}{ value of parameter.}
}

\details{
The distribution has density 
\deqn{f(x) = (x + 1)/(\beta * (\beta + 1)) * exp(-x/\beta).}
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
f <- function(x){dlindley(x, 3)}
## Plot using graphics package
require(graphics)
curve(f, 0, 6)
mode <- molindley(3)
lines(c(mode, mode), c(0, f(mode)), lwd = 2)
}
