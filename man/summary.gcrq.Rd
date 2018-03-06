\name{summary.gcrq}
\alias{summary.gcrq}
%\alias{print.gcrq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarizing model fits for growth charts regression quantiles
}
\description{
summary and print methods for class \code{gcrq}
}
\usage{
\method{summary}{gcrq}(object, digits = max(3, getOption("digits") - 3), ...)

%\method{print}{summary.gcrq}(x,

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class \code{"gcrq"}.
}
  \item{digits}{
controls number of digits printed in output.
}
  \item{\dots}{
further arguments.
}
}
\details{
This method is at a preliminary stage. \code{summary.gcrq} returns some information on the fitted quantile curve at different probability values, such as the estimates, standard errors, values of check (objective) function values at solution. Currently there is no \code{print.summary.gcrq} method, so 
\code{summary.gcrq} itself prints results.
}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Vito M.R. Muggeo
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{gcrq}}
}
\examples{
## see ?gcrq
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ nonlinear }
\keyword{ regression }% __ONLY ONE__ keyword per line
