print.gcrq <-function(x, digits = max(3, getOption("digits") - 4), ...){
#cat("\nNo. of obs:", n, "  No. of estimated parameters:", p, " (",x$pLin+sum(x$info.smooth$dropvcList) ,"unpenalized )" ,"\n")
      n<-nrow(as.matrix(x$fitted.values))
      p<-nrow(as.matrix(x$coefficients))
      n.tau<-length(x$taus)
      cat("Formula: ")
      print(x$call$formula)
      cat("\nQuantile curve at: ", paste(x$taus, collapse="  "), "\n")
      cat("No. of obs:", n, "  No. of param.:", p,"(for each tau);", p*n.tau,"(total)\n")
      if(length(x$info.smooth)>0){
         cat("\nEquivalent Degrees of Freedom:\n")
         if(rownames(x$edf.j)[1]=="Xlin") rownames(x$edf.j)[1]<-"unpenalized"
         attr(x$edf.j,"df.nc")<- NULL
         print(round(x$edf.j, digits))
      }
      all.sic <- log(x$rho/n) + log(n) * colSums(as.matrix(x$edf.j))/(2*n)
      sic<- sum(all.sic) #sum(log(x$rho/n)) +log(n)*sum(x$edf.j)/(2*n)
      cat("\nOverall check =", round(sum(x$rho),digits), "  SIC =", round(sic,digits), " (on edf =",paste(round(sum(x$edf.j),digits-1),")",sep="") ,"\n")
}


