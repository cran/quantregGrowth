summary.gcrq <-
function(object, digits = max(3, getOption("digits") - 3), ...){
      n<-length(attr(object$BB[[1]],"covariate.n"))
      cat("\n*** Noncrossing regression quantiles via P-splines***\n")
      cat("\nCall:\n")
      print(object$call)
      p<-nrow(object$coefficients)
      n.tau<-length(object$taus)
      cat("\nNo. of obs:", n, "  No. of estimated parameters:", p,"(for each curve);", p*n.tau,"(total)\n")
      cat("\nQuantile curves at percentiles: ", object$taus, "\n")
      cat("\nCheck functions:", round(object$rho,2), "\n")
      sic<- sum(log(object$rho/n)) +log(n)*sum(object$df)/(2*n)
      cat("\nOverall Check function =", round(sum(object$rho),digits), "  SIC =", round(sic,digits),"\n")
      }
