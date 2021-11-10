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
\method{summary}{gcrq}(object, type=c("sandw","boot"), digits = max(3, getOption("digits") - 3), 
    signif.stars =getOption("show.signif.stars"), ...)
%\method{print}{summary.gcrq}(x,

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class \code{"gcrq"}.
}
  \item{type}{
Which covariance matrix should be used to compute the estimate standard errors? \code{'boot'} means case-resampling bootstrap (see \code{n.boot} in \code{gcrq()}), \code{'sandw'} mean via the sandwich formula.
}

  \item{digits}{
controls number of digits printed in output.
}
  \item{signif.stars}{
Should significance stars be printed? 
}
  \item{\dots}{
further arguments.
}
}
\details{
\code{summary.gcrq} returns some information on the fitted quantile curve at different probability values, such as the estimates, standard errors, values of check (objective) function values at solution. Currently there is no \code{print.summary.gcrq} method, so 
\code{summary.gcrq} itself prints results.

The SIC returned by \code{print.gcrq} and \code{summary.gcrq} is computed as \eqn{\log(\rho_\tau/n) + \log(n)  edf/(2 n)}{log(rho/n) + log(n) * edf/(2 * n)}, where \eqn{\rho_tau}{rho} is the usual asymmetric sum of residuals (in absolute value). For multiple \eqn{J}{J} quantiles it is \eqn{\log(\sum_\tau\rho_\tau/(n J)) + \log(n J)  edf/(2 n J)}{log(Sum rho/(n*J)) + log(n*J) * edf/(2 * n *J)}. Note that computation of SIC in \code{\link{AIC.gcrq}} relies on the Laplace assumption for the response. 
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
##summary(o)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ nonlinear }
\keyword{ regression }% __ONLY ONE__ keyword per line
