
summary.heckit5 <- function( object) {
   RSq <- function(model, intercept) {
      ## Calculate r-squared.  Note that the way lm() finds R2 is a bit naive -- it checks for intercept
      ## in the formula, but not whether the intercept is present in any of the data vectors (or matrices)
       if(class(model) == "lm") {
          y <- model.response(model.frame(model))
          if(intercept) {
             R2 <- 1 - sum(residuals(model)^2)/sum((y - mean(y))^2)
             R2adj <- 1 - (1 - R2)*(NObs(model) - 1)/(NObs(model) - NParam(model))
          }
          else {
             R2 <- 1 - sum(residuals(model)^2)/sum(y^2)
             R2adj <- 1 - (1 - R2)*(NObs(model))/(NObs(model) - NParam(model))
          }
       }
       else {
          R2 <- model$eq[[ 1 ]]$r2
          R2adj <- model$eq[[ 1 ]]$adjr2
       }
       c(R2, R2adj)
   }
   r <- RSq(object$twoStep1, object$param$oIntercept1)
   R21 <- r[1]
   R2adj1 <- r[2]
   r <- RSq(object$twoStep2, object$param$oIntercept2)
   R22 <- r[1]
   R2adj2 <- r[2]
   iBetaS <- object$param$index$betaS
   stdd <- sqrt(diag(vcov(object, part="full")))
   estimate <- coefTable(coef(object, part="full"), stdd, object$param$df)
   s <- list(estimate=estimate,
             rSquared=list(R21=R21, R2adj1=R2adj1, R22=R22, R2adj2=R2adj2),
             param=object$param)
   class(s) <- c("summary.heckit5", class(s))
   s
}

print.summary.heckit5 <- function( x,
                                 digits=max(3, getOption("digits") - 3),
                                 signif.stars=getOption("show.signif.stars"),
                                  part="full",
                                 ...) {
   cat("2-step estimation of switching regression (tobit-5) model\n")
   cat(x$param$NObs, " observations (", x$param$N1, " selection 1 and ",
       x$param$N2, " selection 2) and ",
       x$param$NParam, " free parameters (df = ",
       x$param$df, ")\n", sep="")
   if(part == "full") {
      i <- x$param$index$betaS
      cat("Probit selection equation:\n")
      printCoefmat(x$estimate[i,], signif.legend=FALSE)
   }
                                        #
   i <- x$param$index$betaO1
   cat("Outcome equation 1:\n")
   printCoefmat(x$estimate[i,], signif.legend=FALSE)
   cat("Multiple R-Squared:", round(x$rSquared$R21, digits),
       ",\tAdjusted R-Squared:", round(x$rSquared$R2adj1, digits), "\n", sep="")
                                        #
   i <- x$param$index$betaO2
   cat("Outcome equation 2:\n")
   if(part == "full")
       printCoefmat(x$estimate[i,], signif.legend=FALSE)
   else
       printCoefmat(x$estimate[i,], signif.legend=TRUE)
   cat("Multiple R-Squared:", round(x$rSquared$R22, digits),
       ",\tAdjusted R-Squared:", round(x$rSquared$R2adj2, digits), "\n", sep="")
                                        #
   if(part=="full") {
      i <- c(x$param$index$sigma1, x$param$index$sigma2, x$param$index$rho1, x$param$index$rho2)
      cat("Error terms:\n")
      printCoefmat(x$estimate[i,], signif.legend=TRUE)
                                        # Here we have a problem -- signif legend is only printed, if
                                        # any of the symbols are printed.  However, this part of the
                                        # vcov is unimplemented and hence it will not print the legend.
   }
   invisible( x )
}
