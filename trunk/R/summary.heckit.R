
summary.heckit <- function( object, ...) {
   s <- list()  # list for results that will be returned

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
   r <- RSq( object$lm, object$param$oIntercept )
   R2 <- r[1]
   R2adj <- r[2]

   stdd <- sqrt(diag(vcov(object, part="full")))
   s$estimate <- coefTable(coef(object, part="full"), stdd, object$param$df)
   s$rSquared <- list(R2=R2, R2adj=R2adj)
   s$param    <- object$param
   class(s) <- c("summary.heckit", class(s))
   return( s )
}

print.summary.heckit <- function( x,
                                 digits=max(3, getOption("digits") - 3),
                                 signif.stars=getOption("show.signif.stars"),
                                 part="full",
                                 ...) {
   cat("2-step estimation of heckit (tobit-2) model\n")
   cat(x$param$NObs, " observations (", x$param$N0, " censored and ", x$param$N1, " observed) and ",
       x$param$NParam, " parameters (df = ",
       x$param$df, ")\n", sep="")
   if(part == "full") {
      i <- x$param$index$betaS
      cat("Probit selection equation:\n")
      printCoefmat(x$estimate[i,], signif.legend=FALSE)
   }
                                        #
   i <- x$param$index$betaO
   cat("Outcome equation:\n")
   if(part == "full")
       printCoefmat(x$estimate[i,], signif.legend=FALSE)
   else
       printCoefmat(x$estimate[i,], signif.legend=TRUE)
   cat("Multiple R-Squared:", round(x$rSquared$R2, digits),
       ",\tAdjusted R-Squared:", round(x$rSquared$R2adj, digits), "\n", sep="")
                                        #
   if(part=="full") {
      i <- c(x$param$index$Mills, x$param$index$sigma, x$param$index$rho)
      cat("Error terms:\n")
      printCoefmat(x$estimate[i,], signif.legend=TRUE)
                                        # Here we have a problem -- signif legend is only printed, if
                                        # any of the symbols are printed.  However, this part of the
                                        # vcov is unimplemented and hence it will not print the legend.
   }
   invisible( x )
}
