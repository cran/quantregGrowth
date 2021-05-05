\name{ncross.rq.fitXB}
\alias{ncross.rq.fitXB}
\alias{ncross.rq.fitX}
%\alias{ncross.rq.fitB}
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
ncross.rq.fitXB(y, x, B=NULL, X=NULL, taus, monotone=FALSE, concave=FALSE, 
    nomiBy=NULL, byVariabili=NULL, ndx=10, deg=3, dif=3, lambda=0, eps=.0001, 
    var.pen=NULL, penMatrix=NULL, lambda.ridge=0, dropcList=FALSE, decomList=FALSE, 
    vcList=FALSE, dropvcList=FALSE, centerList=FALSE, colmeansB=NULL, ...)

ncross.rq.fitX(y, X = NULL, taus, lambda.ridge = 0, eps = 1e-04, 
    sgn.constr=1, adjX.constr=TRUE, ...)

gcrq.rq.cv(y, B, X, taus, monotone, concave, ndx, lambda, deg, dif, var.pen=NULL, 
    penMatrix=NULL, lambda.ridge=0, dropcList=FALSE, decomList=FALSE, 
    vcList=vcList, dropvcList=FALSE, nfolds=10, foldid=NULL, eps=.0001, ...)
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
%  \item{interc}{
%should the model intercept be estimated? see the same argument in \code{\link{gcrq}}
%}
  \item{monotone}{
numerical value (-1/0/+1) to define a non-increasing, unconstrained, and non-decreasing flexible fit, respectively. 
}
  \item{concave}{
numerical value (-1/0/+1) to possibly define concave or convex fits. 
}
  \item{nomiBy}{
useful for VC models (when \code{B} is not provided).
}
  \item{byVariabili}{
useful for VC models (when \code{B} is not provided).
}
  \item{ndx}{
number of internal intervals within the covariate range, see \code{ndx} in \code{\link{ps}}.
}
  \item{deg}{
spline degree, see \code{\link{ps}}.
}
  \item{dif}{
difference order of the spline coefficients in the penalty term.
}
  \item{lambda}{
smoothing parameter value(s), see \code{lambda} in \code{\link{ps}}.
}
  \item{eps}{
tolerance value.
}
  \item{var.pen}{
Varying penalty, see \code{\link{ps}}.
}
  \item{penMatrix}{
Specified penalty matrix, see \code{pen.matrix} in \code{\link{ps}}.
}
%  \item{err.rho}{
% The smoothing parameter is estimated via the ratio \code{err.rho}/\code{u.rho} see \code{\link{gcrq}}.
%}
%  \item{u.rho}{
% The smoothing parameter is estimated via the ratio \code{err.rho}/\code{u.rho} see \code{\link{gcrq}}.
%}
  \item{lambda.ridge}{
a (typically very small) value, see \code{lambda.ridge} \code{\link{gcrq}}.
}
  \item{dropcList}{
see \code{dropc} in \code{\link{ps}}.
}
  \item{decomList}{
see \code{decompose} in \code{\link{ps}}.
}
  \item{vcList}{
to indicate if the smooth is VC or not,  see \code{by} in \code{\link{ps}}.
}
  \item{dropvcList}{
see \code{\link{ps}}.
}
  \item{centerList}{
see \code{center} in \code{\link{ps}}.
}
  \item{colmeansB}{
see \code{center} in \code{\link{ps}}.
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
  \item{sgn.constr}{optional. Vector of signs for the noncrossing constraints. Appropriate only with linear models.}
  \item{adjX.constr}{logical to shift the linear covariates. Appropriate only with linear models.
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
