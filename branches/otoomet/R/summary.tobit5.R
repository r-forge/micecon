summary.tobit5 <- function(object, ...) {
    opt <- summary(object$results)
                                        # class maximisation
    a <- list(opt=opt)
    class(a) <- "summary.tobit5"
    a
}

