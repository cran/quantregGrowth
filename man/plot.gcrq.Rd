\name{plot.gcrq}
\alias{plot.gcrq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot method for gcrq objects
}
\description{
Displaying the estimated growth charts from a \code{gcrq} fit.
}
\usage{
\method{plot}{gcrq}(x, term=NULL, add = FALSE, res = FALSE, conf.level=0, axis.tau=FALSE, 
   interc=TRUE, se.interc=FALSE, legend = FALSE, select.tau, deriv = FALSE, cv = FALSE, 
   transf=NULL, lambda0=FALSE, shade=FALSE, overlap=NULL, rug=FALSE, overall.eff=TRUE, %%n.points=100,
   grid=NULL, smoos=NULL, split=FALSE, shift=0, type=c("sandw","boot"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a fitted \code{"gcrq"} object.
}
  \item{term}{
the variable name or its index in the formula entering the model. Can be vector. Both linear ad spline terms (i.e. included in the model via \code{ps}) can be specified and relevant fitted quantile curves (as optionally specified by \code{select.tau}) will be plotted. If the model includes both linear and smooth terms, the smooth terms are counted and drawn *first*: therefore if the model formula is \code{y~z+ps(x)}, \code{term=1} refers to the smooth term. If \code{NULL}, all smooth terms are plotted according to the \code{split} argument.
If the model includes multiple quantile curves and \code{axis.tau=TRUE}, \code{term=1} refers to the model intercept (if in the model). The variable name should be reported within 'ps()', e.g. \code{'ps(age)'} regardless of additional arguments
specified in \code{ps}.
}
  \item{interc}{
Should the smooth term be plotted along with the model intercept (provided it is included in the model)? Of course such argument is ignored if the smooth term has been called via \code{ps(, dropc=FALSE)} and the plot always includes implicitly the `intercept'. Note that \code{interc=TRUE} is requested to display the noncrossing curves (if multiple quantile curves are being plotted).  
}
  \item{se.interc}{
logical. If \code{TRUE} the standard errors of fitted quantile curves account for uncertainty of the model intercept (provided it is included in the model). If \code{FALSE}, then the uncertainty relates purely to the (usually centred) smooth itself. Ignored if \code{conf.level=0}.
}

  \item{add}{
logical. If \code{TRUE} the fitted quantile curves are added on the current plot.
}
  \item{res}{
logical. If \code{TRUE} `partial residuals' are also displayed on the plot. Borrowing terminology from GLM, partial residuals for covariate 
\eqn{X_j}{Xj} are defined as \cr fitted values corresponding to \eqn{X_j}{Xj} + residuals (from the actual fit). \cr If there is a single covariate, the partial residuals correspond to observed data. If multiple quantile curves have been estimated, the fitted values coming from the `middle' quantile curve are employed to compute the partial residuals. `Middle' means `corresponding to the \eqn{\tau_k}{tau} closest to 0.50'. I don't know if that is the best choice.
}
  \item{conf.level}{
logical. If larger than zero, pointwise confidence intervals for the fitted quantile curve are also shown (at the confidence level specified by \code{conf.level}). Such confidence intervals are independent of the possible intercept accounted for via the \code{intercept} argument. See \code{type}
to select different methods (bootstrap or sandwich) to compute the standard errors. 
}
  \item{axis.tau}{
logical. If \code{TRUE}, the estimated coefficient \code{term} is plotted against the probability values. This graph could be useful if the model has been estimated at several \code{tau} values. 
}
  \item{legend}{
logical. If \code{TRUE} a legend is drawn on on the right side of the plot.
}
  \item{select.tau}{
an optional numeric vector to draw only some of the fitted quantiles. Percentile values or integers 1 to \code{length(tau)} may be supplied. 
}                                                                           
  \item{deriv}{
logical. If \code{TRUE} the first derivative of the fitted curves are displayed.
}
  \item{cv}{
logical. If \code{TRUE} and the \code{"gcrq"} object contains a single smooth term wherein lambda has been selected via CV, then %the matrix \code{cv}, then 
the cross-validation scores against the lambda values are plotted.
}
  \item{transf}{
An optional character string (with \code{"y"} as argument) meaning a function to apply to the predicted values (and possibly residuals) 
  before plotting. E.g. \code{"(exp(y)-0.1)"}. If \code{NULL} (default) it is taken as the inverse of function \code{transf} (*if*) supplied in \code{gcrq}. See argument \code{"transf"} in \code{gcrq()}. If \code{transf} has been specified in \code{gcrq()}, use \code{transf="y"} to force plotting on the transformed scale, i.e. without back transforming.
}
  \item{lambda0}{
logical. If \code{cv=TRUE}, should the CV plot include also the first CV value? Usually the first CV value is at lambda=0, 
and typically it is much bigger than the other values making the plot not easy to read. Default to \code{FALSE} not to display the
first CV value in the plot.
}
  \item{shade}{
logical. If \code{TRUE} and \code{conf.level}>0, the pointwise confidence intervals are portrayed via shaded areas.
}
  \item{overlap}{
NULL or numeric (scalar or vector). If provided and different from \code{NULL}, it represents the abscissa values (on the covariate scale) where the legends (i.e. the probability values) of each curve are set. It will be recycled, if its length differs from the number of quantile curves. If unspecified (i.e. \code{overlap=NULL}), the legends are placed outside the fitted lines on the right side. If specified, \code{legend=TRUE} is implicitly assumed.
}
  \item{rug}{
logical. If \code{TRUE}, the covariate distribution is displayed as a rug plot at the foot of the plot. Default to \code{FALSE}.
}
%  \item{n.points}{
%numeric. Number of values used to plot the fitted curves. Large values provide smoother curves.
%}
%  \item{edf.ylab}{
%Should the edf value to be reported as y label? If \code{NULL}, \code{edf.ylab} is set to TRUE only if there is a single quantile curve to be plotted.
%}
  \item{overall.eff}{
logical. If the smooth term has been called via \code{ps(.., decom=TRUE)}, by specifying \code{overall.eff=TRUE} the overall smooth effect is drawn, otherwise only the penalized part is portrayed (always \emph{without} intercept).
}
  \item{grid}{
if provided, a grid of horizontal and vertical lines is drawn. \code{grid} has to be a list with the following components \code{x,y,col,lty,lwd}. 
If \code{x} (\code{y}) is a vector, the vertical (horizontal) lines are drawn at these locations. If \code{x} (\code{y}) is a scalar, the vertical (horizontal) lines are drawn at \code{x} (\code{y}) equispaced values. \code{col, lty,lwd} refer to the lines to be drawn.
}
  \item{smoos}{
logical, indicating if the residuals (provided that \code{res=TRUE}) will be drawn using a \emph{smoothed} scatterplot. If \code{NULL} (default) the smoothed scatterplot will be employed when the number of observation is larger than 10000.
}
  \item{split}{
logical. If there are multiple terms (both smooth and linear) and \code{split=TRUE}, \code{plot.gcrq()} tries to split the plotting area in 2 columns and number of rows depending on the number of smooths. If \code{split=FALSE}, the plots are produced on the current device according to the current graphics settings. Ignored if there is single smooth term. 
}
  \item{shift}{
Numerical value(s) to be added to the curve(s) to be plotted. If vector with length equal to the number of quantile curves to plot, the \code{shift[j]} is added to the jth quantile curve. 
%Default is \code{NULL} which means the 'model intercept' for VC smooth term and 0 otherwise.
%It can be useful to plot the VC term accounting for the model intercept. 
}
  \item{type}{
If \code{conf.level>0}, which covariance matrix should be used to compute and to portray the pointwise confidence intervals? \code{'boot'} means case-resampling bootstrap (see \code{n.boot} in \code{gcrq()}, \code{'sandw'} mean via the sandwich formula.
}

  \item{\dots}{
Additional graphical parameters:\cr 
\code{xlab}, \code{ylab}, \code{ylim}, and \code{xlim} (effective when \code{add=FALSE});\cr
\code{lwd}, \code{lty}, and \code{col} for the fitted quantile lines; \code{col<0} means color palette for the different curves; \cr
\code{cex} and \code{text.col} for the legend (if \code{legend=TRUE} or \code{overlap} is specified); \cr
\code{cex.p}, \code{col.p}, and \code{pch.p} for the points (if \code{res=TRUE}).\cr
When \code{axis.tau=TRUE}, all arguments accepted by \code{plot()}, \code{points()}, \code{matplot()}, and \code{matpoints()} but \code{pch, type, xlab, ylab, lty}.
}
}
\details{
Takes a "gcrq" object and diplays the fitted quantile curves as a function of the covariate specified in \code{term}. If \code{conf.level}>0 pointwise confidence intervals are also displayed. When the object contains the component \code{cv}, \code{plot.gcrq} can display cross-validation scores against the lambda values, see argument \code{cv}. If a single quantile curve is being displayed, the default 'ylab' includes the relevant edf value (leaving out the basis intercept). If \code{axis.tau=TRUE} and the fit includes several quantile curves, \code{plot.gcrq()} portrays the estimated coefficients versus the probability values. If \code{term} refers to a categorical variable, the point estimates against the categories are plotted (\code{conf.level} is ignored).
% per gli IC, col e' quello dei quantili, lwd = lwd(quantili)*.5, e lty=lty(quantili)+1 (se un solo valore e' specificato) oppure e' proprio quello  
%dei quantili se questo e' un vettore.
}
\value{
The function simply generates a new plot or adds fitted curves to an existing one. 
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Vito M. R. Muggeo %%with contribution by Gianluca Sottile (gianluca.sottile@unipa.it)
}
\note{
Plotting non-crossing curves could depend on the arguments 'interc' and 'shift', in turn depending on how the model has been specified. Take care about that! %%For instance, see fit \code{o2} in Example 3 in \code{\link{gcrq}}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{gcrq}}, \code{\link{predict.gcrq}}
}
\examples{

\dontrun{
## use the fits from ?gcrq
#The additive model
plot(o, res=TRUE, col=2, conf.level=.9, shade=TRUE, split=TRUE)


par(mfrow=c(2,2))
plot(m5, select.tau=c(.1,.5,.9), overlap=0.6, legend=TRUE)
plot(m5, grid=list(x=8,y=5), lty=1) #a 8 times 5 grid.. 
plot(m7, cv=TRUE) #display CV score versus lambda values
plot(m7, res=TRUE, grid=list(x=5, y=8), col=4) #fitted curves at the best lambda value

}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
\keyword{ nonlinear }
