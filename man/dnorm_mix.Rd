\name{Density Normal Mixture}

\alias{dnorm_mix}

\title{ Normal Mixture Density}

\description{ Gives the density value of a univariate Normal Mixture.}

\usage{dnorm_mix(x, pi, mean, sd)}

\arguments{
  \item{x}{ value or vector of values, to apply in the mixture probability density function.}
  \item{pi}{ vector of mixture weights.}
  \item{mean}{ vector of mixture means.}
  \item{sd}{ vector of mixture standard desviations.}
}

\details{
Each component of the mixture has density
\deqn{f_j(x) = 1/(\sqrt(2 * \pi) * \sigma_j) * exp(-(x - \mu_j)/(2 * \sigma_j^(2))}
where \eqn{\mu_j} is the mean of the jth component and \eqn{\sigma_j} is the standard deviation of the jth component. j = 1, 2, ..., g.
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
f <- function(x){dnorm_mix(x, pi = c(0.5, 0.5), mean = c(4, 11), sd = c(2, 2))}
## Plot using graphics package
require(graphics)
curve(f, 0, 15)

require(ggplot2)
## Plot using ggplot2 package
ggplot(data.frame(x=c(0, 15)), aes(x)) +
  stat_function(fun=f, size = 1L, colour = "#0c4c8a") +
  theme_minimal()

## Integrate
integrate(f, "-Inf", "Inf")$value
}
