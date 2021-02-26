\name{vcov.gcrq}
\alias{vcov.gcrq}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Variance-Covariance Matrix for a Fitted 'cgrq' Model
}

\description{
Returns the (currently only bootstrap-based) variance-covariance matrix of the parameter estimates of a fitted gcrq model object.
}

\usage{
\method{vcov}{gcrq}(object, term, type=c("boot","sandw"), ...)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{object}{
a fitted model object of class \code{"gcrq"} returned by \code{gcrq()}.
  }
  \item{term}{
if specified the returned covariance matrix includes entries relevant to parameter estimates for that 'term' only. If missing, the returned matrices refer to all model parameter estimates.
  }
  \item{type}{
Which cov matrix should be returned? Currentl,y only \code{'boot'} is allowed.
  }
  \item{...}{
additional arguments.  
}
}
\details{
If the \code{"gcrq"} object includes results from bootstrap runs (namely the component \code{boot.coef} is not \code{NULL}),
\code{vcov.gcrq()} computes the covariance matrix for the parameter estimates of each quantile curve. If not, 100 bootstrap replicates are run.
}
\value{
A list (with length equal to the length of \code{tau} specified in \code{gcrq}) of square matrices. Namely the list includes the covariance matrices of the parameter estimates for each regression quantile curve.
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
