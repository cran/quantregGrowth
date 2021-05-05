\name{vcov.gcrq}
\alias{vcov.gcrq}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Variance-Covariance Matrix for a Fitted 'gcrq' Model
}

\description{
Returns the variance-covariance matrix of the parameter estimates of a fitted gcrq model object.
}

\usage{
\method{vcov}{gcrq}(object, term, type=c("sandw","boot"), ...)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{object}{
a fitted model object of class \code{"gcrq"} returned by \code{gcrq()}.
  }
  \item{term}{
if specified, the returned covariance matrix includes entries relevant to parameter estimates for that 'term' only. If missing, the returned matrices refer to all model parameter estimates. Currently \code{term} is not allowed. 
  }
  \item{type}{
Which cov matrix should be returned? \code{'boot'} means case-resampling bootstrap (see \code{n.boot} in \code{gcrq()}), \code{'sandw'} mean via the 
sandwich formula.
  }
  \item{...}{
additional arguments.  
}
}
\details{
Bootstrap-based covariance matrix, i.e. \code{type="boot"}, is computable only if the object fit has been obtained by specifying \code{n.boot>0} in \code{gcrq()}.
}
\value{
A list (its length equal the length of \code{tau} specified in \code{gcrq}) of square matrices. Namely the list includes the covariance matrices of the parameter estimates for each regression quantile curve.
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
