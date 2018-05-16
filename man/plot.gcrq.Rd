\name{plot.gcrq}
\alias{plot.gcrq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot method for gcrq objects
}
\description{
Displaying the estimated growth charts from a \code{gcrq} fit.
}
\usage{
\method{plot}{gcrq}(x, term, add = FALSE, res = FALSE, conf.level=0, legend = FALSE, 
    select.tau, deriv = FALSE, cv = FALSE, transf=NULL, lambda0=FALSE, 
    shade=FALSE, overlap=FALSE, rug=FALSE, n.points=100, grid=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a fitted \code{"gcrq"} object.
}
  \item{term}{
the smooth variable name entering the model via \code{ps}. Relevant fitted quantile curves will be plotted.
It may be missing if the model includes a single smooth term.}

  \item{add}{
logical. If \code{TRUE} the fitted quantile curves are added on the current plot.
}
  \item{res}{
logical. If \code{TRUE} `partial' residuals' are also displayed on the plot. Borrowing terminology from GLM, partial residuals for covariate 
\eqn{X_j}{Xj} are defined as \cr fitted values corresponding to \eqn{X_j}{Xj} + residuals. \cr If there is a single covariate, the partial residuals correspond to observed data. If multiple quantile curves have been estimated, the fitted values coming from the `middle' quantile curve are employed to compute the partial residuals. `Middle' means `corresponding to the \eqn{\tau_k}{tau} closest to 0.50'. I don't know if that is the best choice.
}
  \item{conf.level}{
logical. If larger than zero, pointwise confidence intervals for the fitted quantile curve are also shown (at the confidence level specified by \code{conf.level}).
}
  \item{legend}{
logical. If \code{TRUE} a legend is drawn on on the right side of the plot.
}
  \item{select.tau}{
an optional numeric vector to draw only some of the fitted quantiles. Percentile values or integers 1 to \code{length(tau)} may be supplied. 
}                                                                           
  \item{deriv}{
logical. If \code{TRUE} the first derivative of the curve is displayed.
}
  \item{cv}{
logical. If \code{TRUE} and the \code{"gcrq"} object contains the matrix \code{cv}, then 
the cross-validation scores against the lambda values are plotted.
}
  \item{transf}{
An optional character string (with "y" as argument) meaning a function to apply to the response variable 
  before plotting. E.g. \code{"(exp(y)-0.1)"}. See argument \code{"transf"} in \code{gcrq()}.
}
  \item{lambda0}{
logical. If \code{cv=TRUE}, should the CV plot include also the first CV value? Usually the first CV value is at lambda=0, 
and typically it is much bigger than the other values making the plot not easy to read. Default to \code{FALSE} not to display the
first CV value in the plot.
}
  \item{shade}{
logical. If \code{TRUE} and \code{conf.level}>0, the pointwise confidence intervals are portrayed via shaded areas.
}
  \item{overlap}{
logical. If \code{TRUE}, and \code{legend=TRUE}, the legend crosses the fitted curve on the right side. Otherwise the leged is placed outside the 
fitted lines.
}
  \item{rug}{
logical. If \code{TRUE}, the covariate distribution is displayed as a rug plot at the foot of the plot. Default to \code{FALSE}.
}
  \item{n.points}{
numeric. Number of values used to plot the fitted curves. Large values provide smoother curves.
}
  \item{grid}{
if provided, a grid of horizontal and vertical lines is drawn. \code{grid} has to be a list with the following components \code{x,y,col,lty,lwd}. 
If \code{x} (\code{y}) is a vector, the vertical (horizontal) lines are drawn at these locations. If \code{x} (\code{y}) is a scalar, the vertical (horizontal) lines are drawn at \code{x} (\code{y}) equspaced values. \code{col, lty,lwd} refer to the lines to be drawn.
}

  \item{\dots}{
Additional graphical parameters:\cr 
\code{xlab}, \code{ylab}, \code{ylim}, and \code{xlim} (effective when \code{add=FALSE});\cr
\code{lwd}, \code{lty}, and \code{col} for the fitted quantile lines; \code{col} means color palette for the different curves \cr
\code{cex} for the legend (if \code{legend=TRUE}); \cr
\code{cex.p}, \code{col.p}, and \code{pch.p} for the points (if \code{res=TRUE}).
}
}
\details{
Takes a "gcrq" object and diplays the fitted quantile curves. If \code{conf.level}>0 pointwise confidence intervals are also displayed. When the object contains the component \code{cv}, \code{plot.gcrq} can display cross-validation scores against the lambda values, see argument \code{cv}.
% per gli IC, col e' quello dei quantili, lwd = lwd(quantili)*.5, e lty=lty(quantili)+1 (se un solo valore e' specificato) oppure e' proprio quello  
%dei quantili se questo e' un vettore.
}
\value{
The function simply generates a new plot or adds fitted curves to an existing one. 
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Vito M. R. Muggeo with contribution by Gianluca Sottile (gianluca.sottile@unipa.it)
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{gcrq}}
}
\examples{
## see ?gcrq
##plot(m5, conf.level=.8, lty=1)
##plot(m5, lty=1:3, select.tau=c(.1,.5,.9))
##plot(m5, grid=list(x=8,y=5)) #a 8 times 5 grid.. 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
\keyword{ nonlinear }
