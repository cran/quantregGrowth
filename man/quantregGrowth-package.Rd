\name{quantregGrowth-package}
\alias{quantregGrowth-package}
\alias{quantregGrowth}
\docType{package}
\title{
Non-Crossing Additive Regression Quantiles and Non-Parametric Growth Charts. 
}
\description{
Fits non-crossing regression quantiles as a function of linear covariates 
  and smooth terms via P-splines with difference penalties. Random intercepts and selection of linear variables are allowed via the lasso penalties. Estimation of (possibly adaptive) smoothing/tuning parameters (for the spline terms, the random intercepts and the variables selector) are carried out efficientely as part of model fitting.}
\details{
\tabular{ll}{
Package: \tab quantregGrowth\cr
Type: \tab Package\cr
Version: \tab 1.7-0\cr
Date: \tab 2023-07-05\cr
License: \tab GPL\cr
}
Package quantregGrowth allows estimation of growth charts via quantile regression. Given a set of percentiles (i.e. probability values), \code{gcrq} estimates non-crossing quantile curves as a flexible function of quantitative covariates (typically age in growth charts), and possibly additional linear terms. To ensure flexibility, B-splines with a difference \eqn{L_1}{L1} penalty are employed to estimate non parametrically the curves wherein monotonicity and concavity constraints may be also set. Multiple smooth terms, including varying coefficients, are allowed and the amount of smoothness for each term is efficiently included in the model fitting algorithm, see Muggeo et al. (2021). \code{plot.gcrq} displays the fitted lines along with observations and poitwise confidence intervals. 
}
\author{
Vito M.R. Muggeo

Maintainer: Vito M.R. Muggeo <vito.muggeo@unipa.it>
}
\references{
Muggeo VMR, Torretta F, Eilers PHC, Sciandra M, Attanasio M (2021).
Multiple smoothing parameters selection in additive regression quantiles,
\emph{Statistical Modelling}, \bold{21}, 428-448. 

Muggeo VMR (2021). Additive Quantile regression with automatic smoothness selection: the R package quantregGrowth. \cr
https://www.researchgate.net/publication/350844895

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
