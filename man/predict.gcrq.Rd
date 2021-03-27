\name{predict.gcrq}
\alias{predict.gcrq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prediction for "gcrq" objects
}
\description{
Takes a "gcrq" objects and computes fitted values
}

\usage{
\method{predict}{gcrq}(object, newdata, se.fit=FALSE, transf=NULL, xreg, 
                      type=c("sandw","boot"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{object}{
a fitted \code{"gcrq"} object.
}
\item{newdata}{
a dataframe including \emph{all} the covariates of the model. The smooth term is represented by a covariate 
and proper basis functions will be build accordingly. If omitted, the fitted values are used. Ignored if \code{xreg} is provided. 
}
\item{se.fit}{
logical. If \code{TRUE}, standard errors of the fitted quantiles are computed using the bootstrap covariance matrix.
}
\item{transf}{
An optional character string (with \code{"y"} as argument) meaning a function to apply to the predicted values. E.g. \code{"(exp(y)-0.1)"}. 
If \code{NULL} (default) it is taken as the inverse of function \code{transf} (*if*) supplied in \code{gcrq}. The standard errors (provided \code{se.fit=TRUE} has been set) are adjusted accordingly via the Delta method. See argument \code{"transf"} in \code{gcrq()}. If \code{transf} has been specified in \code{gcrq()}, use \code{transf="y"} to force predictions on the transformed scale, i.e. without back transforming.
}
\item{xreg}{
the design matrix for which predictions are requested. Note \code{xreg} has to include the basis functions of the 
B-spline.
}
  \item{type}{
If \code{se.fit=TRUE}, which cov matrix should be used? \code{'boot'} means case-resampling bootstrap (see \code{n.boot} in \code{gcrq()}), \code{'sandw'} mean via the sandwich formula.
}
\item{\dots}{ arguments passed to other functions }
}
\details{
\code{predict.gcrq} computes fitted quantiles as a function of observations included in \code{newdata} or \code{xreg}. 
Either \code{newdata} or \code{xreg} have to be supplied, but \code{newdata} is ignored 
  when \code{xreg} is provided.
}
\value{
If \code{se.fit=FALSE}, a matrix of fitted values with number of rows equal to number of rows of input data 
and number of columns depending on the number of fitted quantile curves (i.e length of \code{taus}). If \code{se.fit=TRUE}, a list of matrices (fitted values and standard errors).
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Vito M.R. Muggeo
}
%\note{
%This function is at a preliminary stage and it should be replaced by the method \code{predict.gcrq}. 
%Please use it with care. This function is not expect to work if the fitted model includes factors or terms like \code{poly()}.
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{gcrq}}, \code{\link{plot.gcrq}}
}
\examples{
##see ?gcrq
## predict(m1, newdata=data.frame(x=c(.3,.7)))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
\keyword{ nonlinear }% __ONLY ONE__ keyword per line
