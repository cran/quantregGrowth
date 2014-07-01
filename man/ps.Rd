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
ps(x, monotone = 0, lambda = 0, pdiff = 3, ndx = NULL, deg = 3, 
    var.pen = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The quantitative covariate supposed to have a nonlinear relationships with the quantiles. 
In growth charts this variable is typically the age.
}
  \item{monotone}{
Numeric value to set up monotonicity restrictions on the fitted smooth function
\itemize{
 \item  '0' = no constrain;
 \item  '1' = non decreasing smooth function;
 \item '-1' = non increasing smooth function.
}
}
  \item{lambda}{
A supplied smoothing parameter for the smooth term. If it is a vector, cross validation is performed to select the `best' 
value. 
}
  \item{pdiff}{
The difference order of the penalty. Default to 3.
}
  \item{ndx}{
The number of intervals of the covariate range used to build the B-spline basis. If \code{NULL}, default, 
the empirical rule of Ruppert is used, namely \eqn{min(n/4,40)}{min(n/4,40)}.
}
  \item{deg}{
The degree of the spline polynomial. Default to 3. 
}
  \item{var.pen}{
A character indicating the varying penalty. See Details.
}
}
\details{

When \code{lambda}=0 an unpenalized fit is obtained. The fit gets smoother as lambda increases, and for a very large value of lambda 
it approaches to a polynomial of degree \code{pdiff-1}. It is also possible to put a varying penalty to set 
a different amount of smoothing. For instance for a 
constant smoothing (\code{var.pen=NULL}) the penalty is \eqn{\lambda\sum_k \Delta_k^2}{lambda sum_k \Delta^2_k} where 
\eqn{\Delta_k}{Delta_k} is the k-th difference (of order \code{pdiff}) of the spline coefficients. When a varying penalty is set, 
the penalty becomes \eqn{\lambda\sum_k \Delta_k^2 w_k}{lambda sum_k \Delta^2_k w_k}. The weights \eqn{w_k}{w_k} depend on \code{var.pen}; 
for instance \code{var.pen="((1:k)^2)"} results in \eqn{w_k=k^2}{w_k=k}. See model \code{m5} in examples of \code{\link{gcrq}}.
}
\value{
The function simply returns the covariate with added attributes relevant to 
smooth term.
}
\references{
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
\keyword{ smooth }% __ONLY ONE__ keyword per line
