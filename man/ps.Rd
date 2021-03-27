\name{ps}
\alias{ps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Specifying a smooth term in the gcrq formula.
}
\description{
Function used to define the smooth term (via P-splines) within the gcrq formula. 
The function actually does not evaluate a (spline) smooth, but simply it 
passes relevant information to proper fitter functions.
     }
\usage{
ps(..., lambda = -1, d = 3, by=NULL, ndx = NULL, deg = 3, knots=NULL,
    monotone = 0, concave = 0, var.pen = NULL, pen.matrix=NULL, dropc=TRUE, 
    center=TRUE, K=2, ridge=FALSE, decompose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{...}{
The quantitative covariate supposed to have a nonlinear relationships with the quantiles. 
In growth charts this variable is typically the age.
}
  \item{lambda}{
A supplied smoothing parameter for the smooth term. If it is negative scalar, the smoothing parameter is estimated iteratively as discussed in Muggeo et al. (2020). If a positive scalar, it represents the actual smoothing parameter. If it is a vector, cross validation is performed to select the `best' value. See Details in \code{\link{gcrq}}. 
}
  \item{d}{
The difference order of the penalty. Default to 3.
}
  \item{by}{
if different from \code{NULL}, a numeric or factor variable of the same dimension as the covariate in \code{...} In the numeric vector case the elements multiply the smooth, evaluated at the corresponding covariate values (i.e. a varying coefficient model).
%Uo' ESSERE un fattore?...........
}
  \item{ndx}{
The number of intervals of the covariate range used to build the B-spline basis. Non-integer values are rounded by \code{round()}. If \code{NULL}, default, the empirical rule of Ruppert is used, namely \eqn{min(n/4,40)}{min(n/4,40)}. It could be reduced (but no less than 10, say) if the sample size is not large and the default value leads to some error in the fitting procedure, see section \code{Note} in \code{\link{gcrq}}.
}
  \item{deg}{
The degree of the spline polynomial. Default to 3. The B-spline basis is composed by \code{ndx}+\code{deg} basis functions and if \code{dropc=TRUE} the first column is removed for identifiability (and the model intercept is estimated without any penalty). 
%However if a varying coefficient term is specified via \code{by}, removing a column is not requested for the identifiability issue. 
}
  \item{knots}{
The knots locations. If \code{NULL}, equispaced knots are set.
}
 \item{monotone}{
Numeric value to set up monotonicity restrictions on the first derivative of fitted smooth function
\itemize{
 \item  '0' = no constraint (default);
 \item  '1' = non-decreasing smooth function;
 \item '-1' = non-increasing smooth function.
}
}
 \item{concave}{
Numeric value to set up monotonicity restrictions on the second derivative of fitted smooth function
\itemize{
 \item  '0' = no constraint (default);
 \item  '1' = concave smooth function;
 \item '-1' = convex smooth function.
}
}
  \item{var.pen}{
A character indicating the varying penalty. See Details.
}
  \item{pen.matrix}{
if provided, a penalty matrix \eqn{A}{A}, say, such that \eqn{A^T A}{A'A} is the penalty matrix actually used in the estimation process. It overwrites \code{ridge}.
}
  \item{dropc}{
logical. Should the first column of the B-spline basis be dropped for the basis identifiability? Default to \code{TRUE}. Note, if \code{dropc=FALSE} is set, 
it is necessary to omit the model intercept AND \code{center=FALSE} or to set \code{lambda.ridge>0} to fit the model.
}
  \item{center}{
logical. If \code{TRUE} the smooth effects are 'centered' over the covariate values, i.e. \eqn{\sum_i \hat{f}(x_i)=0}{sum_i f(xi)=0}.
}
  \item{K}{
A factor tuning selection of wiggliness of the smoothed curve. The larger \code{K}, the smoother the curve. Simulations suggest \code{K=2}. See details.
}
  \item{ridge}{
logical. If \code{TRUE} a ridge penalty is empoyed.
}
  \item{decompose}{
logical. If \code{TRUE}, the B-spline is decomposed into truncated power functions 
such as [x, ..., x^{d-1}, Z]. See Details.
}
}
\details{

When \code{lambda}=0 an unpenalized fit is obtained. At 'middle' lambda values, the fitted curve is a \emph{piecewise} polynomial of order \code{d-1}.
The fit gets smoother as lambda increases, and for a very large value of lambda, it approaches to a polynomial of degree \code{d-1}. 
It is also possible to put a varying penalty to set a different amount of smoothing. Namely for a 
constant smoothing (\code{var.pen=NULL}) the penalty is \eqn{\lambda\sum_k |\Delta^d_k|}{lambda sum_k |\Delta^d_k|} where 
\eqn{\Delta^d_k}{Delta^d_k} is the k-th difference (of order \code{d}) of the spline coefficients. For instance if \eqn{d=1}{d=1}, 
 \eqn{|\Delta^1_k|=|b_k-b_{k-1}|}{|\Delta^1_k|=|b_k-b_{k-1}|} where the \eqn{b_k}{b_k}  are the spline coefficients.
When a varying penalty is set, the penalty becomes \eqn{\lambda\sum_k |\Delta_k^d| w_k}{lambda sum_k |\Delta^d_k| w_k}. 
The weights \eqn{w_k}{w_k} depend on \code{var.pen}; for instance \code{var.pen="((1:k)^2)"} results in \eqn{w_k=k^2}{w_k=k^2}. See model \code{m5} in examples of \code{\link{gcrq}}.
}
\value{
The function simply returns the covariate with added attributes relevant to 
smooth term.
}
\references{
Muggeo VMR, Torretta F, Eilers PHC, Sciandra M, Attanasio M (2020).
Multiple smoothing parameters selection in additive regression quantiles,
\emph{Statistical Modelling}, to appear.


For a general discussion on using B-spline and penalties in regression model see

Eilers PHC, Marx BD. (1996) Flexible smoothing with B-splines and penalties. 
Statistical Sciences, 11:89-121.
}
\author{
Vito M. R. Muggeo
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{gcrq}}
}
\examples{
##see ?gcrq

##gcrq(y ~ ps(x),..) #it works (default: center = TRUE, dropc = TRUE)
##gcrq(y ~ 0 + ps(x, center = TRUE, dropc = FALSE)) #it does NOT work
##gcrq(y ~ 0 + ps(x, center = FALSE, dropc = FALSE)) #it works

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
\keyword{ smooth }% __ONLY ONE__ keyword per line
