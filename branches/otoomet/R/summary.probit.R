summary.probit <- function(object, ...) {
   ## summary for probit -- adds Likelihood Ratio Test to summary.maxLik
   summaryML <- summary.maxLik(object, ...)
   pchi2 <- pchisq(object$LRT$LRT, object$LRT$df, lower.tail=FALSE)
   a <- c(summaryML,
          LRT=list(c(object$LRT, pchi2=pchi2)),
          NParam=object$NParam,
          NObs=object$NObs,
          df=object$df)
   class(a) <- c("summary.probit", class(summaryML))
   a
}
