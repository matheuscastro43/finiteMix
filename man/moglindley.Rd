\name{Mode Generalized Lindley Distribution}

\alias{moglindley}

\title{ Gemeralized Lindley Distribution Mode}

\description{
Gives the mode value of a univariate Generalized Lindley Distribution.
}

\usage{moglindley(alpha, beta, gamma)}

\arguments{
  \item{alpha}{ value of parameter.}
  \item{beta}{ value of parameter.}
  \item{gamma}{ value of parameter.}
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
f <- function(x){dglindley(x, 2, 2, 3)}
## Plot using graphics package
require(graphics)
curve(f, 0, 9)
mode <- moglindley(2, 2, 3)
lines(c(mode, mode), c(0, f(mode)), lwd = 2)
}
