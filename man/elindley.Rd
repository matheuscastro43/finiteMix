\name{Inferences Lindley Distribution}

\alias{elindley}

\title{ Lindley Distribution Estimation}

\description{ Estimates parameters of a univariate Lindley Distribution with maximum likelihood.}

\usage{
elindley(data, plot.it = TRUE, empirical = FALSE, col.estimated = "orange",
         col.empirical = "navy", ...)
}

\arguments{
\item{data}{ vector containing the sample, or list obtained with \link[finiteMix]{rlindley}.}
\item{plot.it}{ logical, TRUE to plot the histogram with estimated distribution curve.}
\item{empirical}{ logical, TRUE to add the empirical curve ("Kernel Density Estimation") in the plot.}
\item{col.estimated}{ a colour to be used in the curve of estimated density.}
\item{col.empirical}{ a colour to be used in the curve of empirical density.}
\item{...}{ further arguments and graphical parameters passed to \link[graphics]{hist}.}
}


%\details{
%%  ~~ If necessary, more details than the description above ~~
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

% \seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{
## Generate a sample.
data = rlindley(n = 1000, beta = 2)

## And now, estimate the parameters, using the 'data' list.
elindley(data)

## Not plotting the graphic.
elindley(data, plot.it = FALSE)

## Adding the empirical curve to the graphic.
elindley(data, empirical = TRUE)

## Changing the color of the curves.
elindley(data, empirical = TRUE, col.estimated = "pink", 
         col.empirical = "red3")

## Using "..."
elindley(data, empirical = TRUE, col.estimated = "pink", 
         col.empirical = "red3", breaks = 300)}