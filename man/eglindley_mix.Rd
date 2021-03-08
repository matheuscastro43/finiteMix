\name{Inferences Generalized Lindley Mixture}

\alias{eglindley_mix}

\title{ Generalized Lindley Mixture Estimation}

\description{ Estimates parameters of a univariate Generalized Lindley Mixture with k-means clustering and EM-algorithm.}

\usage{
eglindley_mix(data, g, lim.em = 100, criteria = "dif.psi", plot.it = TRUE,
              empirical = FALSE, col.estimated = "orange", col.empirical = "navy",
              ...)
}

\arguments{
\item{data}{ vector containing the sample, or list obtained with \link[finiteMix]{rglindley_mix}.}
\item{g}{ number of components in the mixture.}
\item{lim.em}{ limit of EM Iterations.}
\item{criteria}{ the stop criteria to be used, could be "dif.psi" to calculate differences on
parameters matrix or "dif.lh" to calculate differences on Likelihood function.}
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
data = rglindley_mix(n = 1000, pi = c(0.6, 0.4), alpha = c(10, 18), beta = c(1, 2), gamma = c(3, 4))

## And now, estimate the parameters, using the 'data' list.
eglindley_mix(data, g = 2, lim.em = 5)

## Or using the sample vector.
# eglindley_mix(data$sample, g = 2, lim.em = 5)

## Using the diference in the log-likelihood as stop criteria.
# eglindley_mix(data, g = 2, criteria = "dif.lh", lim.em = 10)

## Not plotting the graphic.
# eglindley_mix(data, g = 2, plot.it = FALSE, lim.em = 5)

## Adding the empirical curve to the graphic.
#eglindley_mix(data, g = 2, empirical = TRUE, lim.em = 5)

## Changing the color of the curves.
# eglindley_mix(data, g = 2, lim.em = 5, empirical = TRUE, col.estimated = "pink",
#               col.empirical = "red3")

## Using "..."
# eglindley_mix(data, g = 2, lim.em = 5, empirical = TRUE, col.estimated = "pink",
#               col.empirical = "red3", breaks = 300)}