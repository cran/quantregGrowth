## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  ##collapse = TRUE,
  comment = ">",
  fig.align="center",
  out.extra='style="display:block; margin: auto"',
  tidy=FALSE
)

## ----1, message = F, fig.width=5, fig.height=4--------------------------------
library(quantregGrowth)
data(SiChildren)

o <-gcrq(height ~ ps(age), tau=seq(.05,.95,l=7), data=SiChildren)

plot(o, res=TRUE) #res=TRUE displays data too

## -----------------------------------------------------------------------------
charts(o, k=c(10,10.5,11,16,17)) #the quantile at the specified k values

## ---- echo=FALSE, results='asis'----------------------------------------------
knitr::kable(charts(o, k=19:21))

## ----1a, message = F----------------------------------------------------------
oM<-gcrq(height ~ ps(age, monotone=1), data=SiChildren, tau=seq(.05,.95,l=7))

## ----1b, message = F, fig.width=7, fig.height=4-------------------------------
par(mfrow=c(1,2))
#the 1st plot..
plot(o, xlim=c(10.2,12.5), col=1, lty=3, lwd=2)
plot(oM, add=TRUE, col=2, lty=1, lwd=2)
#the 2nd..
plot(oM, legend=TRUE, overlap=15, grid=list(x=15,y=10), col=2, lty=1, 
     ylab="Height (cm)", xlab="Age (years)")

## ----2, message=FALSE---------------------------------------------------------
set.seed(1515)
d<-mgcv::gamSim(n=200, eg=1, verbose=FALSE) #verbose=TRUE suppresses the message..

o <- gcrq(y ~ ps(x0) + ps(x1) + ps(x2) + ps(x3), data=d, tau=.5)

## ----fig-margin, fig.width=7, fig.height=5------------------------------------
# Plot the fit
par(mfrow=c(2,2))
plot(o, res=TRUE, col=2, conf.level=.95, shade=TRUE, cex.p=.6) #cex.p<1 to reduce the points..


## ----3, message = F-----------------------------------------------------------
o1<-gcrq(y ~ ps(x0) + x1 + ps(x2) + x3, data=d, tau=.5)

#the summary method
summary(o1)


## ----4, message = F-----------------------------------------------------------
o2<-gcrq(y ~ x1 + x3 + ps(x0) + ps(x2), data=d, tau=c(.12,.25,.5,.7,.8,.9))


## -----------------------------------------------------------------------------
n=100
x<-1:n/n
y<-2*x+rnorm(n)*rev(x)

o3 <-gcrq(y ~ x, tau=seq(.05,.95,l=10)) #fit 10 (noncrossing) regression quantiles

## ---- fig.width=7, fig.height=3-----------------------------------------------
par(mfrow=c(1,3))
plot(x,y)
plot(o3, add=TRUE)
plot(o3, term=1, axis.tau=TRUE)
plot(o3, term="x", axis.tau=TRUE, conf.level = .95, col=2)

## -----------------------------------------------------------------------------
data(growthData)

D1 <- diff(diag(23),dif=1)[,-1] #the 1st order diff matrix
w <- c(rep(1,15),rep(1000,7)) #the weight vector s.t. length(w)=nrow(D1)
A <- diag(w) %*% D1 #the penalty matrix

o4 <-gcrq(y~ps(x, ndx=20, pen.matrix=A), data=growthData, tau=.5)


## ---- fig.width=7, fig.height=3-----------------------------------------------
o5 <-gcrq(y~ps(x, d=1), data=growthData, tau=.5)
o6 <-gcrq(y~ps(x, d=1, monotone = 1), data=growthData, tau=.5)

par(mfrow=c(1,3))
plot(o4, res=TRUE, col=2, lwd=3, ylim=c(0,12))
plot(o5, res=TRUE, col=3, lwd=3, ylim=c(0,12))
plot(o6, res=TRUE, col=4, lwd=3, ylim=c(0,12))

## ---- echo=FALSE, results='asis'----------------------------------------------
knitr::kable(head(mtcars, 10))
