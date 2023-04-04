\name{charts}
\alias{charts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Easy computing growth charts
}
\description{
Computes and returns quantiles as a function of the specified covariate values
}
\usage{
charts(fit, k, file = NULL, digits=2, conf.level=0, 
    transf=NULL, se.type=c("sandw","boot"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
The object fit returned by \code{\link{gcrq}}
}
  \item{k}{
Scalar or vector indicating the covariate values. If scalar, \code{k} equispaced values in the covariate range are taken. 
}
  \item{file}{
If specified, the (path) file name wherein the returned matrix including the quantiles will be written via \code{write.csv()}
}
  \item{digits}{
Number of digits whereby the estimated quantiles are rounded. 
}
  \item{conf.level}{
If larger than zero, the pointwise confidence intervals for the estimated quantiles are also returned. If \code{conf.level=0} the simple point estimates. 
}
  \item{transf}{
An optional character string (with \code{"y"} as argument) meaning a function to apply to the predicted values. E.g. \code{"(exp(y)-0.1)"}. 
If \code{NULL} (default) it is taken as the inverse of function \code{transf} (*if*) supplied in \code{gcrq}. The standard errors (provided \code{se.fit=TRUE} has been set) are adjusted accordingly via the Delta method. See argument \code{"transf"} in \code{gcrq()}. If \code{transf} has been specified in \code{gcrq()}, use \code{transf="y"} to force predictions on the transformed scale, i.e. without back transforming.

}
  \item{se.type}{
Which covariance matrix should be used, provided that \code{conf.level>0}. See \code{type} in \code{\link{predict.gcrq}}. 
}
  \item{\dots}{
Further arguments passed on to \code{write.csv()}
}
}
\details{
This function is simply a wrapper for \code{\link{predict.gcrq}}
}
\value{
A matrix having number of columns equal to the number of quantile curves and number of rows depending \code{k} 
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Vito Muggeo
}
\note{
\code{charts} just works with models having a single smooth term. See \code{\link{predict.gcrq}} when the model involves multiple covariates.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{predict.gcrq}}
}
\examples{
\dontrun{
charts(_fit_, k=10)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
