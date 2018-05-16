\name{quantregGrowth-package}
\alias{quantregGrowth-package}
\alias{quantregGrowth}
\docType{package}
\title{
Growth Charts via Regression Quantiles
}
\description{
Fits non-crossing regression quantiles as a function of linear covariates 
  and smooth terms via B-splines with difference penalties.}
\details{
\tabular{ll}{
Package: \tab quantregGrowth\cr
Type: \tab Package\cr
Version: \tab 0.4-1\cr
Date: \tab 2018-05-16\cr
License: \tab GPL\cr
}
Package quantregGrowth allows estimation of growth charts via quantile regression. Given a set of percentiles, \code{gcrq}
estimates non-crossing quantile curves as a flexible function of quantitative covariates (typically age), and possibly 
additional linear terms. To ensure flexibility, B-splines with a difference \eqn{L_1}{L1} penalty are employed to estimate nonparametrically 
the curves;  additionally monotonicity constraints may be also set. \code{plot.gcrq} displays the fitted lines along with observations and poitwise confidence intervals. 
}
\author{
Vito M.R. Muggeo

Maintainer: Vito M.R. Muggeo <vito.muggeo@unipa.it>
}
\references{

Muggeo VMR, Sciandra M, Tomasello A, Calvo S (2013).
Estimating growth charts via nonparametric quantile
regression: a practical framework with application
in ecology, \emph{Environ Ecol Stat}, \bold{20}, 519-531.

Muggeo VMR (2018). Using the R package quantregGrowth: some examples. \cr
https://www.researchgate.net/publication/323573492

Some references on growth charts (the first two papers employ the so-called LMS method)

Cole TJ, Green P (1992) Smoothing reference centile curves: the LMS method and penalized likelihood.
\emph{Statistics in Medicine} \bold{11}, 1305-1319.

Rigby RA, Stasinopoulos DM (2004) Smooth centile curves for skew and kurtotic data modelled using the
Box-Cox power exponential distribution. \emph{Statistics in Medicine} \bold{23}, 3053-3076.

Wei Y, Pere A, Koenker R, He X (2006) Quantile regression methods for reference growth charts. 
\emph{Statistics in Medicine} \bold{25}, 1369-1382.



Some references on regression quantiles

Koenker R (2005) Quantile regression. Cambridge University Press, Cambridge.

Cade BS, Noon BR (2003) A gentle introduction to quantile regression for ecologists. \emph{Front Ecol Environ}
\bold{1}, 412-420.


}

\keyword{ package }
\keyword{ models }
\keyword{ regression }


\seealso{
\code{\link{gcrq}}, \code{\link[quantreg]{rq}} in package quantreg
}
\examples{
#see ?gcrq for some examples
}
