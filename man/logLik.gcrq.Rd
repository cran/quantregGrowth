\name{logLik.gcrq}
\alias{logLik.gcrq}
\alias{AIC.gcrq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Log Likelihood, AIC and BIC for gcrq objects
}
\description{
The function returns the log-likelihood value(s) evaluated at the estimated coefficients
}
\usage{
\method{logLik}{gcrq}(object, summ=TRUE, ...)
\method{AIC}{gcrq}(object, ..., k=2)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
A \code{gcrq} fit returned by \code{gcrq()}
}
  \item{summ}{
If \code{TRUE}, the log likelihood values (and relevant edf) are summed over the different taus to provide a unique value accounting for the different quantile curves. If \code{FALSE}, tau-specific values are returned.
}
\item{k}{
Optional numeric specifying the penalty of the edf in the AIC formula. \code{k < 0} means \code{k=log(n)}.
}
  \item{\dots}{
optional  arguments (nothing in \code{logLik.gcrq}). For \code{AIC.gcrq}, \code{summ=TRUE} or \code{FALSE} can be set. 
}
}
\details{
The 'logLikelihood' is computed by assuming an asymmetric Laplace distribution for the response as in \code{\link{logLik.rq}}, namely \eqn{n (\log(\tau(1-\tau))-1-\log(\rho_\tau/n))}{n * (log(tau * (1 - tau)) - 1 - log(object$rho/n))}, where \eqn{\rho_\tau}{object$rho} is the minimized objective function. When there are multiple quantile curves \eqn{j=1,2,...,J} (and \code{summ=TRUE}) the formula is

\eqn{n (\sum_j\log(\tau_j(1-\tau_j))-J-\log(\sum_j\rho_{\tau_j}/(n J)))}{n * ( Sum_j log(tau_j * (1 - tau_j)) - J - log(Sum_j object$rho_j/(n*J)))}

\code{AIC.gcrq} simply returns \code{-2*logLik + k*edf} where \code{k} is 2 or \code{log(n)}.  
}
\value{
The log likelihood(s) of the model fit \code{object}
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
\code{\link{logLik.rq}}
}
\examples{
   \donttest{
## logLik(o) #a unique value (o is the fit object  from gcrq)
## logLik(o, summ=FALSE) #vector of the log likelihood values
## AIC(o, k=-1) #BIC
   }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression} % use one of  RShowDoc("KEYWORDS")
