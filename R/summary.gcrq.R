summary.gcrq <- function(object, type = c("sandw", "boot"), digits = max(3, getOption("digits") - 3), 
         signif.stars = getOption("show.signif.stars"), ...) {
   edfS <- object$edf.j
   if (length(object$info.smooth) > 0) {
      if (rownames(edfS)[1] == "Xlin") 
         edfS <- edfS[-1, , drop = FALSE]
      # rownames(edfS)<-paste('ps(', names(object$BB), ')',sep='')
      colnames(edfS) <- rep("edf", ncol(edfS))
   }
   # cat('\n***Noncrossing regression quantiles via P-splines***\n')
   cat("Call:\n")
   print(object$call)
   
   n <- nrow(as.matrix(object$fitted.values))
   p <- nrow(as.matrix(object$coefficients))
   n.tau <- ncol(as.matrix(object$coefficients))
   #sic <- sum(log(object$rho/n)) + log(n) * sum(object$edf.j)/(2 * n)
   sic<- log(sum(object$rho)/(n*n.tau))+sum(object$edf.j)*log(n*n.tau)/(2*n*n.tau)
   
   list.vcov <- vcov.gcrq(object, type = type)
   for (j in 1:n.tau) {
      est <- as.matrix(object$coefficients)[, j]
      se <- sqrt(diag(list.vcov[[j]]))
      ris <- cbind(Est = est, StErr = se, 
            `|z|` = round(abs(est)/se, 2), `p-value` = pchisq((est/se)^2, df = 1, lower.tail = FALSE))
      rownames(ris) <- rownames(as.matrix(object$coefficients))
      if (length(object$info.smooth) > 0) {
         nomi.p.spline <- unlist(unname(sapply(object$BB[names(object$Bderiv)], function(.x) attr(.x, "coef.names"))))
         nomi.param <- setdiff(rownames(ris), nomi.p.spline)
         ris <- ris[nomi.param, , drop = FALSE]
         cat("\n--------  Percentile:", sprintf("%.2f",object$taus[j]), "  check function: ", formatC(round(object$rho[j], digits - 1), digits=5), "-----\n")
         if (object$pLin > 0) {
            cat("\nparametric terms:\n")
            if(j<n.tau) printCoefmat(ris, digits = digits, signif.stars = signif.stars, signif.legend=FALSE) else 
               printCoefmat(ris, digits = digits, signif.stars = signif.stars)
         }
      cat("\nsmooth terms:\n")
      printCoefmat(edfS[, j, drop = FALSE], digits = digits)
      } else {
         #cat("\n") 
         cat("\n--------  Percentile:", sprintf("%.2f",object$taus[j]), "  check function: ", formatC(round(object$rho[j], digits - 1), digits=5), "-----\n")
         cat("\n")
         if(j<n.tau) printCoefmat(ris, digits = digits, signif.stars = signif.stars, signif.legend=FALSE) else 
            printCoefmat(ris, digits = digits, signif.stars = signif.stars)
      }
   }
   #cat("\n=========================\n\nNo. of obs:", n, "  Check function =", 
   #    round(sum(object$rho), digits - 1), " SIC =", round(sic, digits), 
   #    " ( on edf =", round(sum(object$edf.j), 3), ")", "\n")
   #cat("No. of params:", p, "(for each curve);", p * n.tau, "(total)\n")
   #cat("\n\n")
   
   #nc<-if(sum(object$DF.POS, object$DF.NEG)>0) TRUE else FALSE
   nc <- attr(object$edf.j, "df.nc") & (length(object$tau)>1)
   cat("=========================\n\nNo. of obs =", n, "   No. of params =", 
       p* n.tau, paste("(", p, sep=""),"for each quantile)", "\n")
   cat("Overall check =", round(sum(object$rho), digits - 1), " SIC =", round(sic, digits-1), 
       "on edf =", round(sum(object$edf.j), 2), "(ncross constr:",paste(nc, ")",sep=""), "\n")   

}