\name{predictQR}
\alias{predictQR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prediction for "gcrq" objects
}
\description{
Takes a "gcrq" objects and computes fitted values
}
\usage{
predictQR(object, newdata, xreg)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{object}{
a fitted \code{"gcrq"} object.
}
\item{newdata}{
a dataframe including \emph{all} the covariates of the model. The smooth term is represented by a covariate 
and proper basis functions will be build accordingly. Ignored if \code{xreg} is provided. 
}
  \item{xreg}{
the design matrix for which predictions are requested. Note \code{xreg} has to include the basis functions of the 
B-spline.
}
}
\details{
\code{predictQR} computes fitted quantiles as a function of observations included in \code{newdata} or \code{xreg}. 
Either \code{newdata} or \code{xreg} have to be supplied, but \code{newdata} is ignored 
  when \code{xreg} is provided.
}
\value{
A matrix of fitted values with number of rows equal to number of rows of input data 
and number of columns depending on the fitted quantile curves.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Vito M.R. Muggeo
}
\note{
This function is at a preliminary stage and it should be replaced by the method \code{predict.gcrq}. 
Please use it with care.
}

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
\keyword{ nonlinear }% __ONLY ONE__ keyword per line
