summary.probit <- function(object, ...) {
    opt <- summary(object$results)
                                        # class maximisation
    pchi2 <- pchisq(object$LRT$LRT, object$LRT$df, lower.tail=FALSE)
    a <- list(opt=opt, LRT=c(object$LRT, pchi2=pchi2))
    class(a) <- "summary.probit"
    a
}
