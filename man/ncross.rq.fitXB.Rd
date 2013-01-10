\name{ncross.rq.fitXB}
\alias{ncross.rq.fitXB}
\alias{ncross.rq.fitX}
\alias{ncross.rq.fitB}
\alias{gcrq.rq.cv}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of noncrossing regression quantiles with monotonicity restrictions.
}
\description{
These are internal functions of package \code{quantregGrowth} and should be not 
called by the user.
}

\usage{
ncross.rq.fitXB(y, x, B = NULL, X = NULL, taus, monotone = FALSE, 
  adj.middle = FALSE, ndx = 10, lambda = 0, deg = 3, dif = 3, eps = 1e-04, 
  plott = 0, var.pen = NULL, ...)

ncross.rq.fitX(y, X = NULL, taus, lambda = 0, adj.middle = FALSE, 
  eps = 1e-04, ...)

ncross.rq.fitB(y, x, B = NULL, taus, monotone = FALSE, adj.middle = 
  FALSE, ndx = 10, lambda = 0, deg = 3, dif = 3, plott = 0, 
  var.pen = NULL, eps = 1e-04, ...)

gcrq.rq.cv(y, B, X, taus, monotone, ndx, lambda, deg, dif, var.pen = NULL, 
    cv = TRUE, nfolds = 10, foldid = NULL)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{the responses vector. see \code{\link{gcrq}}
  }
  \item{x}{
the covariate supposed to have a nonlinear relationship.
}
  \item{B}{
the B-spline basis.
}
  \item{X}{
the design matrix for the linear parameters.
}
  \item{taus}{
the percentiles of interest.
}
  \item{monotone}{
numerical value (-1/0/+1) to define a non-increasing, unconstrained, and non-decreasing flexible fit, respectively. 
}
  \item{adj.middle}{
ignore it!
}
  \item{ndx}{
number of internal intervals within the covariate range, see \code{\link{ps}}.
}
  \item{lambda}{
smoothing parameter, see \code{\link{ps}}.
}
  \item{deg}{
spline degree, see \code{\link{ps}}.
}
  \item{dif}{
difference order of the spline coefficients in the penalty term.
}
  \item{eps}{
tolerance value.
}
  \item{plott}{
Plotting results, unuseful. See \code{\link{plot.gcrq}}.
}
  \item{var.pen}{
Varying penalty, see \code{\link{ps}}.
}
  \item{foldid}{
vector (optional) to perform cross validation, see the same arguments in \code{\link{gcrq}}.
}
  \item{nfolds}{
number of folds for crossvalidation, see the same arguments in \code{\link{gcrq}}.
}
  \item{cv}{
returning cv scores; see the same arguments in \code{\link{gcrq}}.
}
  \item{\dots}{
optional.
}
}
\details{
These functions are called by \code{\link{gcrq}} to fit growth charts based on regression 
quantiles with non-crossing and monotonicity restrictions. The computational methods are based on the package 
quantreg by R. Koenker and details are described in the reference paper.
}
\value{
A list of fit information.
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
##See ?gcrq
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ nonlinear }
\keyword{ regression }% __ONLY ONE__ keyword per line
