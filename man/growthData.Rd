\name{growthData}
\alias{growthData}
\docType{data}
\title{Simulated data to illustrate capabilities of the package}
\description{
  The \code{growthData} data frame has 200 rows and 3 columns.
}
\usage{data(growthData)}
\format{
  A data frame with 200 observations on the following 3 variables.
  \describe{
    \item{\code{x}}{the supposed `age' variable.}
    \item{\code{y}}{the supposed growth variable (e.g. weight).}
    \item{\code{z}}{an additional variable to be considered in the model.}
  }
}
\details{
Simulated data to illustrate capabilities of the package.
}
%\source{
%  Bacon D.W., Watts D.G. (1971) Estimating the transistion between two intersecting straight lines.
%  \emph{Biometrika} \bold{58}: 525 -- 534.
%}
%\references{
%  PhD thesis by R.A. Cook
%}
\examples{
data(growthData)
with(growthData, plot(x,y))
}
\keyword{datasets}
