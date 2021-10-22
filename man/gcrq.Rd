\name{gcrq}
\alias{gcrq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Growth charts regression quantiles with automatic smoothness estimation
}
\description{
Modelling unspecified nonlinear relationships between covariates and quantiles of the response conditional distribution. Typical example is estimation nonparametric growth charts (via quantile regression). Quantile curves are estimated via B-splines with 
a \eqn{L_1}{L1} penalty on the spline coefficient differences, while non-crossing and possible monotonicity and concavity restrictions are set to obtain estimates
more biologically plausible. Linear terms can be specified in the model formula. Multiple smooth terms, including varying coefficients, with automatic selection of corresponding smoothing parameters are allowed.  
}
\usage{
gcrq(formula, tau=c(.1,.25,.5,.75,.9), data, subset, weights, na.action, 
    transf=NULL, y=TRUE, n.boot=0, eps=0.001, display=FALSE, 
    method=c("REML","ML"), df.opt=2, df.nc=FALSE, lambda0=.1, h=0.8, lambda.max=2000, 
    tol=0.01, it.max=20, single.lambda=TRUE, foldid=NULL, nfolds=10, 
    lambda.ridge=0, sgn.constr=NULL, adjX.constr=TRUE, contrasts=NULL, ...)
}
\arguments{
  \item{formula}{
a standard R formula to specify the response in the left hand side, and the covariates in the right hand side, such as \code{y~ps(x)+z}, see Details for further eamples.
}
  \item{tau}{
a numeric vector to specify the quantile curves of interest. Default to probability values \eqn{(.1,.25,.5,.75,.9)}.
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
  before fitting. E.g. if \code{y>=0}, it could be advisable to model \code{"log(y+0.1)"}. It can be useful to guarantee fitted values within a specified range. If provided,
  the resulting object fit refer to the model for the transformed response and it will include the corresponding inversefunction (numerically computed) to be used to back transform predictions (see argument \code{transf} in \code{predict.gcrq} and \code{plot.gcrq}).
}
  \item{y}{
logical. If \code{TRUE} (default) the returned object includes also the responses vector.
}
%  \item{interc}{
%logical. \code{interc=TRUE} (default) means that a model intercept is estimated (and the B-spline bases are modified) to allow for a non singular design matrix.
%If this is the case, a very small ridge-type penalty is exploited to allow estimation with a design matrix containing both 
%a full B-spline basis and a column of ones. \code{interc=TRUE} overwrites the intercept specification in the formula (e.g., \code{~0+..}), 
%and it is ignored if the model does not include a \code{"ps"} term.
%}
  \item{n.boot}{
Number of nonparametric (cases resampling) bootstrap samples to be used. If \code{n.boot>0}, the covariance matrix can be obtained as empirical covariance matrix of the bootstrap distributions, see \code{\link{vcov.gcrq}}. Notice that the smoothing parameter (if relevant) is assumed fixed. Namely it does change throughout the bootstrap replicates.  
}
  \item{eps}{
A small positive constant to ensure noncrossing curves (i.e. the minimum distance between two consecutive curves). Use it at your risk! If \code{eps} is large, the resulting fitted quantile curves could appear unreasonable.
}
  \item{display}{
Logical. Should the iterative process be printed? Ignored if no smooth is specified in the formula or if all the smoothing parameters specified in \code{ps} terms are fixed.
}
  \item{method}{
character, \code{"ML"} or \code{"REML"} affecting the smoothing parameter estimation. Default is \code{"REML"} which appears to provide better performance in simulation studies. 
Ignored if no smoothing parameter has to be estimated. 
}
  \item{df.opt}{
How the model and term-specific degrees of freedom are computed. \code{df.opt=1} means via the null penalized coefficients, and \code{df.opt=2} via the trace of the approximate hat  matrix. Ignored if no smoothing parameter is be estimated.
}
  \item{df.nc}{
logical. If \code{TRUE} and the model refers to multiple quantile curves, the degrees of freedom account for the noncrossing constraints. 
Ignored for single quantile fits. Default to \code{FALSE}, as it is still experimental.
}
  \item{lambda0}{
the starting value for the lambdas to be estimated. Ignored if all the smoothing parameters specified in \code{ps} terms are fixed.
}
  \item{h}{
The step halving factor affecting estimation of the smoothing parameters. Lower values lead to slower updates in the lambda values. Ignored if all the smoothing parameters specified in \code{ps} terms are fixed.
}
  \item{lambda.max}{
The upper bound for lambda estimation. Ignored if all the smoothing parameters specified in \code{ps} terms are fixed.
}
  \item{tol}{
The tolerance value to declare convergence. Ignored if all the smoothing parameters specified in \code{ps} terms are fixed.
}
  \item{it.max}{
The maximum number of iterations in lambdas estimation. Ignored if all the smoothing parameters specified in \code{ps} terms are fixed.
}
  \item{single.lambda}{
Logical. Should the smoothing parameter (for each smooth term) to be the same across the quantile curves being estimated? Ignored when just a single quantile curve is being estimated.  
%Ignored if all the smoothing parameters specified in \code{ps} terms are fixed.
}
  \item{foldid}{
optional. A numeric vector identifying the group labels to perform cross validation to select the smoothing parameter. 
Ignored if the \code{lambda} argument in \code{ps()} is not a vector.
}
  \item{nfolds}{
optional. If \code{foldid} is not provided, it is scalar specifying the number of `folds' (groups) which should be used to perform cross validation to select 
the smoothing parameter. Default to 10, but it is ignored if the \code{lambda} argument in \code{ps()} is not a vector.
}
  \item{lambda.ridge}{
Numerical value (typically very small) to stabilize model estimation.
}
  \item{sgn.constr}{optional. Vector of signs for the noncrossing constraints affecting the slopes of linear covariates. If provided, its lenght should be ugual to the number of linear coefficients otherwise it will be recycled. If \code{NULL}, its value is determined by 
  a preliminary heteroscedastic model. Appropriate only with linear models (i.e. no \code{ps} term in the formula).
}
  \item{adjX.constr}{logical. If \code{TRUE}, each linear covariate is shifted (by adding or substracting its min or max) in order make the constraints on the intercept more effective to prevent crossing of quantile curves. 
  Appropriate only with linear models.
}
  \item{contrasts}{an optional list. See argument \code{contrasts.arg} in \code{model.matrix.default}.}
  \item{\dots}{further arguments.}
}
\details{
The function fits regression quantiles at specified percentiles given in \code{tau} as a function of 
covariates specified in the \code{formula} argument. The \code{formula} may include linear terms and  one or several \code{\link{ps}} 
terms to model nonlinear relationships with quantitative covariates, usually age in growth charts. When the \code{lambda} argument 
in \code{ps()} is a negative scalar, the smoothing parameter is estimated iteratively as discussed in Muggeo et al. (2020). If a positive scalar, it represents the actual smoothing parameter value. 

Smoothing parameter selection via 'K-fold' cross validation (CV) is also allowed (but not recommended) if the model includes a single \code{ps} term:  \code{lambda} should be a vector of candidate values, and the final fit is returned at the `optimal' lambda value. To select the smoothing parameter via CV, \code{foldid} or \code{nfolds} may be supplied. If provided \code{foldid} overwrites \code{nfolds}, otherwise \code{foldid} is obtained via random extraction, namely \code{sample(rep(seq(nfolds), length = n))}. However selection of smoothing parameter via CV is allowed only with a unique \code{ps} term in the formula.
}
\value{
This function returns an object of class \code{gcrq}, that is a list with the following components (only the most important are listed)
  \item{coefficients}{The matrix of estimated regression parameters; the number of columns equals the number of the fitted quantile curves.}
  \item{x}{the design matrix of the final fit (including the dummy rows used by penalty).}
  \item{edf.j}{a matrix reporting the edf values for each term at each quantile curve. See the section 'Warning' below.  }
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
  \item{y}{the response vector (if \code{gcrq()} has been called with \code{y=TRUE}).}
  \item{contrasts}{the contrasts used, when the model contains a factor.}
  \item{xlevels}{the levels of the factors (when included) used in fitting.}
  \item{taus}{a vector of values between 0 and 1 indicating the estimated quantile curves.}
  \item{call}{the matched call.}
}

\section{Warning}{
The options \code{'REML'} or \code{'ML'} of the argument \code{method}, refer to how the degrees of freedom are computed to update the lambda estimates. 

Currently, standard errors are obtained via the sandwich formula or the nonparametric bootstrap (case resampling). Both methods ignore uncertainty in the smoothing parameter selection. 

Since version 1.2-1, computation of the approximate edf can account for the noncrossing constraints by specifying \code{df.nc=TRUE}. That could affect model estimation when the smoothing parameter(s) have to be estimated, because the term specific edf are used to update the lambda value(s). When lambda is not being estimated (it is fixed or there is no \code{ps} term in the formula), parameter estimate is independent of the \code{df.nc} value. The \code{summary.gcrq} method reports if the edf account for the noncrossing constraints.     

Using \code{ps(.., center=TRUE)} in the formula leads to lower uncertainty in the fitted curve while guaranteeing noncrossing constraints.  

Currently, decomposition of Bsplines (i.e. \code{ps(..,decom=TRUE)}) is incompatible with shape (monotonicity and concavity) restrictions and even with noncrossing constraints.

%\cr
}
%The function (and underlying method) works pretty well in obtaining point estimates and displaying quantile curves accordingly.
%Typically this is the main (and unique) goal when dealing with growth charts. However from a statistical viewpoint 
%there are some important limitations affecting the theory and the relevant package,
%\enumerate{
%\item Computation of model degrees of freedom
%\item Computation of standard errors, currently implemented via the nonparametric (cases resampling) bootstrap (if \code{n.boot>0})
%  }
%Currently the function does not return standard errors for the parameter estimates (unless \code{n.boot}>0) and degrees of freedom 
%are roughtly computed by summing the 'zero' residuals. Due to noncrossing constraints, the number of zero residuals 
%might not be equal to the number of estimated parameters even in unpenalized models - except for the median, or for models
%where a single quantile curve has been fitted. 
%}

\note{
This function is based upon the package quantreg by R. Koenker.
Currently methods specific to the class \code{"gcrq"} are \code{print.gcrq}, \code{summary.gcrq}, \code{vcov.gcrq}, \code{plot.gcrq}, and \code{predict.gcrq}.

If the sample is not large, and/or the basis rank is large (i.e. a large number of columns) and/or there are relatively few distinct values in the covariate distribution, the fitting algorithm may fail returning error messages like the following 

\code{> Error info =  20 in stepy2: singular design}

To remedy it, it suffices to change some arguments in \code{ps()}: to decrease \code{ndx} or \code{deg} (even by a small amount) or 
to increase (even by a small amount) the lambda value. Sometimes even by changing slightly the tau probability value (for instance from 0.80 to 0.79) can bypass the aforementioned errors.
}

\author{ Vito M. R. Muggeo, \email{vito.muggeo@unipa.it} }

\references{

V.M.R. Muggeo, F. Torretta, P.H.C. Eilers, M. Sciandra, M. Attanasio (2020).
Multiple smoothing parameters selection in additive regression quantiles,
Statistical Modelling, to appear.

V. M. R. Muggeo (2021). Additive Quantile regression with automatic smoothness selection: the R package quantregGrowth. \cr
https://www.researchgate.net/publication/350844895

V. M. R. Muggeo, M. Sciandra, A. Tomasello, S. Calvo (2013).
Estimating growth charts via nonparametric quantile
regression: a practical framework with application
in ecology, Environ Ecol Stat, 20, 519-531.

V. M. R. Muggeo (2018). Using the R package quantregGrowth: some examples. \cr
https://www.researchgate.net/publication/323573492

}


\seealso{
 \code{\link{ps}}, \code{\link{plot.gcrq}}, \code{\link{predict.gcrq}}
}
\examples{
\dontrun{
#An additive examples.. from ?mgcv::gam
d<-mgcv::gamSim(n=200, eg=1)
o<-gcrq(y ~ ps(x0) + ps(x1)+ ps(x2) + ps(x3), data=d, tau=.5, n.boot=50)
plot(o, res=TRUE, col=2, conf.level=.9, shade=TRUE, split=TRUE)


#some simple examples involving just a single smooth
data(growthData) #load data
tauss<-seq(.1,.9,by=.1) #fix the percentiles of interest

m1<-gcrq(y~ps(x), tau=tauss, data=growthData) #lambda estimated..

m2<-gcrq(y~ps(x, lambda=0), tau=tauss, data=growthData) #unpenalized.. very wiggly curves
#strongly penalized models
m3<-gcrq(y~ps(x, lambda=1000, d=2), tau=tauss, data=growthData) #linear 
m4<-gcrq(y~ps(x, lambda=1000, d=3), tau=tauss, data=growthData) #quadratic  

#penalized model with monotonicity restrictions
m5<-gcrq(y~ps(x, monotone=1, lambda=10), tau=tauss, data=growthData)

#monotonicity constraints,lambda estimated, and varying penalty
m6<-gcrq(y~ps(x, monotone=1, lambda=10, var.pen="(1:k)"), tau=tauss, data=growthData) 
m6a<-gcrq(y~ps(x, monotone=1, lambda=10, var.pen="(1:k)^2"), tau=tauss, data=growthData) 

par(mfrow=c(2,3))
plot(m1, pch=20, res=TRUE)
plot(m2, pch=20, res=TRUE)
plot(m3, add=TRUE, lwd=2)
plot(m4, pch=20, res=TRUE)
plot(m5, pch=20, res=TRUE, legend=TRUE, col=2)
plot(m6, lwd=2, col=3)
plot(m6a, lwd=2, col=4)

#select lambda via 'K-fold' CV (only with a single smooth term)
m7<-gcrq(y~ps(x, lambda=seq(0.02,50,l=20)), tau=tauss, data=growthData) 
par(mfrow=c(1,2))
plot(m7, cv=TRUE) #display CV score versus lambda values
plot(m7, res=TRUE, grid=list(x=5, y=8), col=4) #fit at the best lambda (by CV) 


#===VC example

n=50
x<-1:n/n
y0<-10+sin(2*pi*x)
y1<-seq(7,11,l=n)
y<-c(y0,y1)+rnorm(2*n)*.1
x<-c(x,x)
z<-rep(0:1, each=n)

#include as a continuous covariate..
o<-gcrq(y~ ps(x, ndx=5) + ps(x,by=z,ndx=5), n.boot=0, tau=.5)
predict(o, newdata=data.frame(x=c(.3,.7), z=c(0,1)))

#.. as a factor
g<-factor(z)
o1<-gcrq(y~ 0+ps(x,by=g,ndx=5), tau=.5) #you have to remove the model intercept!
predict(o1, newdata=data.frame(x=c(.3,.7), g=factor(c(0,1))))


}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
\keyword{ model }% __ONLY ONE__ keyword per line
\keyword{ smooth }% __ONLY ONE__ keyword per line
