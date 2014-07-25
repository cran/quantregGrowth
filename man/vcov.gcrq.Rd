\name{vcov.gcrq}
\alias{vcov.gcrq}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Variance-Covariance Matrix for a Fitted 'cgrq' Model
}

\description{
Returns the variance-covariance matrix of the parameter estimates of a fitted gcrq model object.
}

\usage{
\method{vcov}{gcrq}(object, ...)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{object}{
a fitted model object of class \code{"gcrq"} returned by \code{gcrq()}.
  }
  \item{...}{
additional arguments.  
}
}
\details{
If the \code{"gcrq"} object includes results from bootstrap runs (namely the component \code{boot.coef} is not \code{NULL}),
\code{vcov.gcrq()} computes the covariance matrix for teh parameter estimates of each quantile curve .
}
\value{
A list including the covariance matrices of the parameter estimates for each regression quantile curve.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Vito Muggeo
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{summary.gcrq}}
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
