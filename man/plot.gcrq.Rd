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
\method{plot}{gcrq}(x, add = FALSE, y = FALSE, legend = FALSE, select.tau, 
    deriv = FALSE, cv = FALSE, transf=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a fitted \code{"gcrq"} object.
}
  \item{add}{
logical. If \code{TRUE} the fitted quantile curves are added on the current plot.
}
  \item{y}{
logical. If \code{TRUE} raw data are also displayed, provided that the object has been called with the argument \code{y=TRUE}.
}
  \item{legend}{
logical. If \code{TRUE} a legend is drawn on on the right side of the fitted curves.
}
  \item{select.tau}{
an optional numeric vector to draw only some of the fitted quantiles. Percentile values or integers 1 to length(tau) may be supplied. 
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
  \item{\dots}{
Additional graphical parameters, such as `xlab', `ylab', and `xlim'; 'lwd', 'col' and 'lty' for the fitted quantile lines; `cex' for 
the legend (if \code{legend=TRUE}); 'cex.p', 'col.p', and 'pch.p' for the points (if \code{y=TRUE}).
}
}
\details{
Takes a "gcrq" object and diplays the fitted quantile curves. When the object contains the component \code{cv}, 
\code{plot.gcrq} can display cross-validation scores against the lambda values, see argument \code{cv}.
}
\value{
The function simply generates a new plot or adds fitted curves to an existing one. 
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Vito M. R. Muggeo
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
\keyword{ nonlinear }
