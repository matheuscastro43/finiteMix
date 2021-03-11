\name{Number of Components}

\alias{g.search}

\title{ Choose the Number of Components on a Finite Mixture}

\description{
Search for the best number of components with non-parametric and parametric criterions
}

\usage{
g.search(data, g = NULL, family = NULL, lim.em = 40, plot.it = TRUE, 
         col.wss = "navy", col.ic = c("red", "orange"), ...)}

\arguments{
  \item{data}{ vector or list containing the sample.}
  \item{g}{ vector with the numbers of components.}
  \item{family}{ parametric family to be used.}
  \item{lim.em}{ limit of EM Iterations.}
  \item{plot.it}{ logical; if TRUE (default) a graph is shown.}
  \item{col.wss}{ a colour to be used in the wss curve.}
  \item{col.ic}{ a colour to be used in the information criterion curve.}
  \item{...}{ further arguments and graphical parameters passed to \link[base]{plot}.}
}

%\details{
%}

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
## Generating a sample
data = rnorm_mix(2000, c(.5, .5), c(1, 9), c(1, 2))
## The histogram of the sample
data$plot

## Finding the number of components
g.search(data, family = "Normal")
}
