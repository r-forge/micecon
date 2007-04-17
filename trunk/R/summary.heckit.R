
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
   if( object$tobitType == 2 ) {
      r <- RSq( object$lm, object$param$oIntercept )
      R2 <- r[1]
      R2adj <- r[2]
      s$rSquared <- list(R2=R2, R2adj=R2adj)
   } else {
      r <- RSq(object$lm1, object$param$oIntercept1)
      R21 <- r[1]
      R2adj1 <- r[2]
      r <- RSq(object$lm2, object$param$oIntercept2)
      R22 <- r[1]
      R2adj2 <- r[2]
      iBetaS <- object$param$index$betaS
      s$rSquared <- list(R21=R21, R2adj1=R2adj1, R22=R22, R2adj2=R2adj2)
   }

   stdd <- sqrt(diag(vcov(object, part="full")))
   s$estimate <- coefTable(coef(object, part="full"), stdd, object$param$df)
   s$param    <- object$param
   s$tobitType <- object$tobitType
   class(s) <- "summary.heckit"
   return( s )
}

print.summary.heckit <- function( x,
                                 digits=max(3, getOption("digits") - 3),
                                 part="full",
                                 ...) {

   cat("--------------------------------------------\n")
   cat("Tobit", x$tobitType, "model" )
   if( x$tobitType == 2 ) {
      cat( " (sample selection model)\n" )
   } else {
      cat( " (switching regression model)\n" )
   }
   cat( "2-step Heckman / heckit estimation\n" )

   cat( x$param$NObs, "observations" )
   if( x$tobitType == 2 ) {
      cat( " (", x$param$N0, " censored and ", x$param$N1, " observed)",
         sep = "" )
   } else {
      cat( " (", x$param$N1, " selection 1 and ",
         x$param$N2, " selection 2)", sep = "" )
   }
   cat( " and", x$param$NParam, "free parameters" )
   cat( " (df = ", x$param$df, ")\n", sep="")
   if(part == "full") {
      i <- x$param$index$betaS
      cat("Probit selection equation:\n")
      printCoefmat(x$estimate[i,], signif.legend=FALSE)
   }
   if( x$tobitType == 2 ) {
      cat("Outcome equation:\n")
      printCoefmat( x$estimate[ x$param$index$betaO, ],
         signif.legend = ( part != "full" ) )
      cat("Multiple R-Squared:", round(x$rSquared$R2, digits),
         ",\tAdjusted R-Squared:", round(x$rSquared$R2adj, digits), "\n", sep="")
   } else {
      cat("Outcome equation 1:\n")
      printCoefmat( x$estimate[ x$param$index$betaO1, ],
         signif.legend = FALSE )
      cat("Multiple R-Squared:", round(x$rSquared$R21, digits),
          ",\tAdjusted R-Squared:", round(x$rSquared$R2adj1, digits), "\n", sep="")
      cat("Outcome equation 2:\n")
      printCoefmat( x$estimate[ x$param$index$betaO2, ],
         signif.legend = ( part != "full" ) )
      cat("Multiple R-Squared:", round(x$rSquared$R22, digits),
         ",\tAdjusted R-Squared:", round(x$rSquared$R2adj2, digits), "\n", sep="")
   }
   if(part=="full") {
      if( x$tobitType == 2 ) {
         i <- c(x$param$index$Mills, x$param$index$sigma, x$param$index$rho)
      } else {
         i <- c( x$param$index$Mills1, x$param$index$Mills2,
            x$param$index$sigma1, x$param$index$sigma2,
            x$param$index$rho1, x$param$index$rho2 )
      }
      cat("Error terms:\n")
      printCoefmat(x$estimate[i,])
                                        # Here we have a problem -- signif legend is only printed, if
                                        # any of the symbols are printed.  However, this part of the
                                        # vcov is unimplemented and hence it will not print the legend.
   }
   cat("--------------------------------------------\n")
   invisible( x )
}
