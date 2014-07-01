summary.gcrq <-
function(object, digits = max(3, getOption("digits") - 3), ...){
      n<-length(attr(object$BB[[1]],"covariate.n"))
      cat("\n***Noncrossing regression quantiles via P-splines***\n")
      cat("\nCall:\n")
      print(object$call)
      n.tau<-length(object$taus)
      if(n.tau>1) {
          p<-nrow(object$coefficients)
          cat("\nNo. of obs:", n, "  No. of estimated parameters:", p,"(for each curve);", p*n.tau,"(total)\n")
          cat("Quantile curves (",length(object$taus),") at percentiles: ", object$taus, "\n")
          cat("Check functions:", round(object$rho,2), "\n")
          sic<- sum(log(object$rho/n)) +log(n)*sum(object$df)/(2*n)
          cat("Overall Check function =", round(sum(object$rho),digits), "  SIC =", round(sic,digits),"\n")
          } else {
          p<-length(object$coefficients)
          cat("\nNo. of obs:", n, "  No. of estimated parameters:", p,"\n")
          cat("Quantile curves at percentile: ", object$taus, "\n")
          sic<- sum(log(object$rho/n)) +log(n)*sum(object$df)/(2*n)
          cat("Check function =", round(sum(object$rho),digits), "  SIC =", round(sic,digits),"\n")
          }
      }
