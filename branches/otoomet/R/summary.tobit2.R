summary.tobit2 <- function(object, ...) {
    opt <- summary(object$results)
                                        # class maximisation
    a <- list(opt=opt)
    class(a) <- "summary.tobit2"
    a
}
