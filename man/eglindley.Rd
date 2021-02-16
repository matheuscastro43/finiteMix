\name{Inferences Generalized Lindley Distribution}

\alias{eglindley}

\title{ Generalized Lindley Distribution Estimation}

\description{ Estimates parameters of a univariate Generalized Lindley Distribution with k-means clustering and EM-algorithm.}

\usage{
eglindley(data, plot.it = TRUE, empirical = FALSE,
                     col.estimated = "orange", col.empirical = "navy", ...)
}

\arguments{
\item{data}{ vector containing the sample, or list obtained with \link[finiteMix]{rglindley}.}
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
data = rglindley(n = 1000, alpha = 18, beta = 2, gamma = 4)

## And now, estimate the parameters, using the 'data' list.
eglindley(data)

## Or using the sample vector.
eglindley(data$sample)

## Not plotting the graphic.
eglindley(data, plot.it = FALSE)

## Adding the empirical curve to the graphic.
eglindley(data, empirical = TRUE)

## Changing the color of the curves.
eglindley(data, empirical = TRUE, col.estimated = "pink", 
          col.empirical = "red3")

## Using "..."
eglindley(data, empirical = TRUE, col.estimated = "pink", 
          col.empirical = "red3", breaks = 300)}