\name{gcrq}
\alias{gcrq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Growth charts regression quantiles
}
\description{
Estimation of nonparametric growth charts via quantile regression. Quantile curves are estimated via B-splines with 
a L1 penalty on the spline coefficient differences, while non-crossing and monotonicity restrictions are set to obtain estimates
more biologically plausible. Linear terms are allowed in the model specification.
}
\usage{
gcrq(formula, tau = c(0.1, 0.25, 0.5, 0.75, 0.9), data, subset, weights, 
      na.action, transf=NULL, y = TRUE, interc=FALSE, 
      foldid = NULL, nfolds = 10, cv = TRUE, n.boot=0, eps=.0001, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
a standard R formula to specify the response in the left hand side, and the covariates in the right hand side. See Details.
}
  \item{tau}{
a numeric vector to specify the percentiles of interest. Default to probability values \eqn{(.1,.25,.5,.75,.9)}.
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
  \item{transf}{
an optional character string (with \code{"y"} as argument) meaning a function to apply to the response variable 
  before fitting. E.g. \code{"log(y+0.1)"}.
}
  \item{y}{
logical. If \code{TRUE} (default) the returned object includes also the responses vector.
}
  \item{interc}{
logical. If \code{formula} includes a \code{"ps"} term, \code{interc=TRUE} means that a model intercept is also estimated.
If this is the case, a very small ridge-type penalty is exploited to allow estimation with a design matrix containing both 
a full B-spline basis and a column of ones. \code{interc=TRUE} overwrites the intercept specification in the formula (e.g., \code{~0+..}), 
and it is ignored if the model does not include a \code{"ps"} term.
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
logical. If \code{TRUE} (default) and there is a single smooth term wherein \code{lambda} is a vector, the returned object includes also the matrix \code{cv} having number of rows equal to length of \code{lambda} and number of columns equal to \code{nfolds}.
Ignored if the \code{lambda} argument in \code{ps()} is \emph{not} a vector. 
}
  \item{n.boot}{
Number of nonparametric (cases resampling) bootstrap samples to be used. Notice that the smoothing parameter (if relevant) does change 
throughout the bootstrap replicates.  Set \code{n.boot>0} if you plan to plot the fitted quantiles along with pointwise confidence intervals.
}
  \item{eps}{
A small positive constant to ensure noncrossing curves. Use it at your risk! If \code{eps} is large, the resulting fitted quantile 
curves could appear unreasonable.
}
  \item{\dots}{
further arguments.
}
}
\details{
The function fits regression quantiles at specified percentiles given in \code{tau} as a function of 
covariates specified in the \code{formula} argument. The \code{formula} may optionally include several \code{\link{ps}} 
terms to model nonlinear relationships with quantitative covariates, usually age in growth charts. When the \code{lambda} argument 
in \code{ps()} is scalar, it represents the actual smoothing 
parameter. When it is a vector, 'K-fold' cross validation is performed to select the `optimal' lambda value and the model is 
fitted at such selected lambda value. To select the smoothing parameter via CV, \code{foldid} or \code{nfolds} 
may be supplied. If provided \code{foldid} overwrites \code{nfolds}, otherwise \code{foldid} is obtained via random extraction, namely 
\code{sample(rep(seq(nfolds), length = n))}. However selection of smoothing parameter is allowed only with a unique ps() term in the formula.
}
\value{
This function returns an object of class \code{gcrq}, that is a list with the following components
  \item{coefficients}{The matrix of estimated regression parameters; the number of columns equals the number of the fitted quantile curves.}
  \item{x}{the design matrix of the final fit (including the dummy rows used by penalty).}
  \item{df}{a vector reporting the df values for each quantile curve. See the section 'Warning' below.  }
  \item{rho}{a vector including the values of the objective functions at the solution for each quantile curve.}

  \item{fitted.values}{a matrix of fitted quantiles (a column for each \code{tau} value)}
  \item{residuals}{a matrix of residuals (a column for each \code{tau} value)}
  \item{D.matrix}{the penalty matrix (multiplied by the smoothing parameter value).}
  \item{D.matrix.nolambda}{the penalty matrix.}
  \item{pLin}{number of linear covariates in the model.}

  \item{info.smooth}{some information on the smoothing term (if included in the formula via \code{ps}).}
  \item{BB}{further information on the smoothing term (if present in the formula via \code{ps}), including 
    stuff useful for plotting via \code{plot.gcrq()}.}
  \item{Bderiv}{if the smooth term is included, the first derivative of the B spline basis.}
  \item{boot.coef}{The array including the estimated coefficients at different bootstrap samples (provided that \code{n.boot}>0 has been set).}
  \item{y}{the response vector (if gcrq() has been called with \code{y=TRUE}).}
  \item{contrasts}{the contrasts used, when the model contains a factor.}
  \item{xlevels}{the levels of the factors (when included) used in fitting.}
  \item{taus}{a vector of values between 0 and 1 indicating the estimated quantile curves.}
  \item{call}{the matched call.}
}
\references{
V. M. R. Muggeo, M. Sciandra, A. Tomasello, S. Calvo (2013).
Estimating growth charts via nonparametric quantile
regression: a practical framework with application
in ecology, Environ Ecol Stat, 20, 519-531.

V. M. R. Muggeo (2018). Using the R package quantregGrowth: some examples. 
https://www.researchgate.net/publication/323573492

}
\author{ Vito M. R. Muggeo, \email{vito.muggeo@unipa.it} }

\note{
This function is based upon the package quantreg by R. Koenker.
Currently methods specific to the class \code{"gcrq"} are \code{print.gcrq}, \code{summary.gcrq}, \code{vcov.gcrq}, and \code{plot.gcrq}.
}

\section{Warning}{
The function (and underlying method) works pretty well in obtaining point estimates and displaying quantile curves accordingly.
Typically this is the main (and unique) goal when dealing with growth charts. However from a statistical viewpoint 
there are some important limitations affecting the theory and the relevant package,
\enumerate{
\item Computation of model degrees of freedom
\item Computation of standard errors, currently implemented via the nonparametric (cases resampling) bootstrap (if \code{n.boot>0})
}
Currently the function does not return standard errors for the parameter estimates (unless \code{n.boot}>0) and degrees of freedom 
are roughtly computed by summing the 'zero' residuals. Due to noncrossing constraints, the number of zero residuals 
might not be equal to the number of estimated parameters even in unpenalized models - except for the median, or for models
where a single quantile curve has been fitted. 
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
m2<-gcrq(y~ps(x, mon=0, lambda=1000, d=2), tau=tauss, data=growthData) #linear 
m3<-gcrq(y~ps(x, mon=0, lambda=1000, d=3), tau=tauss, data=growthData) #quadratic  

#penalized model with monotonicity restrictions
m4<-gcrq(y~ps(x, mon=1, lambda=10), tau=tauss, data=growthData)

#monotonicity constraints with varying penalty
m5<-gcrq(y~ps(x, mon=1, lambda=10, var.pen="(1:k)^3"), tau=tauss, data=growthData) 

par(mfrow=c(2,2))
plot(m1, pch=20, res=TRUE)
plot(m2, pch=20, res=TRUE)
plot(m3, add=TRUE, lwd=2)
plot(m4, pch=20, res=TRUE)
plot(m5, pch=20, res=TRUE, legend=TRUE)

#select lambda via 'K-fold' CV
m6<-gcrq(y~ps(x, lambda=seq(0,100,l=20)), tau=tauss, data=growthData) 
par(mfrow=c(1,2))
plot(m6, cv=TRUE) #display CV score versus lambda values
plot(m6, res=TRUE) #fitting at the best lambda value

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
\keyword{ model }% __ONLY ONE__ keyword per line
\keyword{ smooth }% __ONLY ONE__ keyword per line
