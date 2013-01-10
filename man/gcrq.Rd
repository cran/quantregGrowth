\name{gcrq}
\alias{gcrq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Growth charts regression quantiles
}
\description{
Estimation of nonparametric growth charts via quantile regression. Quantile curves are estimated via B-splines with 
a quadratic penalty on the spline coefficient differences, and non-crossing and monotonicity restrictions are set to obtain estimates
more biologically plausible. Linear terms are allowed in the model specification.
}
\usage{
gcrq(formula, tau = c(0.1, 0.25, 0.5, 0.75, 0.9), data, subset, weights, 
  na.action, y = TRUE, foldid = NULL, nfolds = 10, cv = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
a standard R formula to specify the response in the left hand side, and the covariates in the right hand side. See Details.
}
  \item{tau}{
a numeric vector to specify the percentiles of interest. Default to \eqn{(.1,.25,.5,.75,.9)}.
}
  \item{data}{
the dataframe where the variables required by the formula, subset and weights arguments are stored.
}
  \item{subset}{
optional. A vector specifying a subset of observations to be used in the fitting process.
}
  \item{weights}{
optional. A numeric vector specifying weights to be assigned to the observations in the fitting process.
Currently unimplemented. 
}
  \item{na.action}{
a function which indicates how the possible `NA's are handled.
}
  \item{y}{
if \code{TRUE} (default) the returned object includes also the responses vector.
}
  \item{foldid}{
optional. A numeric vector identifying the group labels to perform cross validation to select the smoothing parameter. 
Ignored if the \code{lambda} argument in \code{ps()} is not a vector.
}
  \item{nfolds}{
optional. If \code{foldid} is not provided, it is scalar specifying the number of `folds' (groups) which should be used to perform cross validation to select 
the smoothing parameter. Default to 10, but it is ignored if the \code{lambda} argument in \code{ps()} is not a vector.
}
  \item{cv}{
logical. If \code{TRUE} and the \code{lambda} argument in \code{ps()} is a vector, the returned object includes also the matrix \code{cv} 
having number of rows equal to length of \code{lambda} and number of columns equal to \code{nfolds}. 
}
  \item{\dots}{
further arguments.
}
}
\details{
The function fits regression quantiles at specified percentiles given in \code{tau} as a function of 
covariates specified in the \code{formula} argument. The \code{formula} may optionally include a \code{\link{ps}} 
term to model a nonlinear relationship with a quantitative covariate, usually age. When the \code{lambda} argument in \code{ps()} is scalar, it represents the actual smoothing 
parameter. When it is a vector, cross validation is performed to select the `optimal' lambda value: at this aim \code{foldid} or \code{nfolds} 
may be supplied. If provided \code{foldid} overwrites \code{nfolds}, otherwise \code{foldid} is obtained via random extraction, namely 
\code{sample(rep(seq(nfolds), length = n))}
}
\value{
This function returns an object of class \code{gcrq}, that is a list with the following components
  \item{coefficients}{The matrix of estimated regression parameters; the number of columns equals the number of the fitted quantile curves.}
  \item{B}{the design matrix of the final fit.}
  \item{df}{a vector reporting the df values for each quantile curve. See the section 'Warning' below.  }
  \item{rho}{a vector including the values of the objective functions at the solution for each quantile curve.}
  \item{BB}{some information on the smoothing term (if included in the formula via \code{ps}).}
  \item{BB}{further information on the smoothing term (if included in the formula via \code{ps}).}
  \item{Bderiv}{if the smooth term is included, the first derivative of the B spline basis.}
  \item{y}{the response vector (if gcrq() has been called with \code{y=TRUE}).}
  \item{contrasts}{the contrasts used, when the model contains a factor.}
  \item{xlevels}{the levels of the factors (when included) used in fitting.}
  \item{taus}{a vector of values between 0 and 1 indicating the quantile curves to be estimated.}
  \item{call}{the matched call.}
}
\references{
V. M. R. Muggeo, M. Sciandra, A. Tomasello, S. Calvo (2013).
Estimating growth charts via nonparametric quantile
regression: a practical framework with application
in ecology, Environ Ecol Stat, to appear.
}
\author{ Vito M. R. Muggeo, \email{vito.muggeo@unipa.it} }

\note{
This function is based upon the package quantreg by R. Koenker.
Currently methods specific to the class \code{"gcrq"} are \code{plot.gcrq} and \code{summary.gcrq}
}

\section{Warning}{
The function (and underlying method) works pretty well in obtaining point estimates and displaying quantile curves accordingly.
Typically this is the main (and unique) goal when dealing with growth charts. However from a statistical viewpoint 
there are some important limitations affecting the theory and the relevant package,
\enumerate{
\item Computation of model degrees of freedom
\item Computation of standard errors
}
Currently the function does not return standard errors for the parameter estimates and degrees of freedom are roughtly computed 
by summing the zero residuals for model containing the smooth term, or simply by the number of parameters in linear models.
}

\seealso{
 \code{\link{ps}}, \code{\link{plot.gcrq}}
}
\examples{
\dontrun{
data(growthData) #load data
tauss<-seq(.1,.9,by=.1) #fix the percentiles of interest

m1<-gcrq(y~ps(x, mon=0), tau=tauss, data=growthData) #unpenalized.. very wiggly curves
#strongly penalized models
m2<-gcrq(y~ps(x, mon=0, lambda=1000, pdiff=2), tau=tauss, data=growthData) #lines 
m3<-gcrq(y~ps(x, mon=0, lambda=1000, pdiff=3), tau=tauss, data=growthData) #quadratic  

#penalized model with monotonicity restrictions
m4<-gcrq(y~ps(x, mon=1, lambda=10), tau=tauss, data=growthData)

#monotonicity constraints with varying penalty
m5<-gcrq(y~ps(x, mon=1, lambda=10, var.pen="(1:k)^3"), tau=tauss, data=growthData) 

par(mfrow=c(2,2))
plot(m1, pch=20, y=TRUE)
plot(m2, pch=20, y=TRUE)
plot(m3, add=TRUE, lwd=2)
plot(m4, pch=20, y=TRUE)
plot(m5, pch=20, y=TRUE, legend=TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
\keyword{ model }% __ONLY ONE__ keyword per line