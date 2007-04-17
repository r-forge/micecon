summary.selection <- function(object, ...) {
   ## object      object of class "selection" (tobit2 and tobit5)
   ## ...         additional arguments for "summary.maxLik"
   ## 
   ## RESULTS:
   ## list of class "summary.selection" with following components:
   ## maximum    : function value at optimum
   ## estimateProbit : probit part of the estimate
   ## estimateEquation : equation -"-
   ## estimateEps      : sigma & rho part
   ## + additional components of "summary.maxLik"
   ##
   if( "tobit2" %in% class( object ) ) {
      model <- 2
   } else if( "tobit5" %in% class( object ) ) {
      model <- 5
   } else {
      stop( "summary.selection currently works for tobit2 and tobit5",
         " models only" )
   }

   sl <- summary.maxLik(object, ...)
   sl$estimateS <- sl$estimate[object$param$index$betaS,]
   sl$NObs <- object$param$NObs
   sl$N1   <- object$param$N1
   sl$NXS  <- object$param$NXS
   sl$df   <- object$param$df

   if( model == 2 ){
      sl$estimateO <- sl$estimate[object$param$index$betaO,]
      sl$estimateErr <- sl$estimate[c(object$param$index$sigma,
                                         object$param$index$rho),]
      sl$N0   <- object$param$N0
      sl$NXO  <- object$param$NXO
   } else if( model == 5 ){
      sl$estimateO1 <- sl$estimate[object$param$index$betaO1,]
      sl$estimateO2 <- sl$estimate[object$param$index$betaO2,]
      sl$estimateErr <- sl$estimate[c(object$param$index$sigma1,
                                         object$param$index$sigma2,
                                         object$param$index$rho1,
                                         object$param$index$rho2),]
      sl$N2   <- object$param$N2
      sl$NXO1 <- object$param$NXO1
      sl$NXO2 <- object$param$NXO2
   }

   class(sl) <- c( "summary.selection", paste( "summary.tobit", model,
      sep = "" ), class(sl) )
   return( sl )
}

print.summary.selection <- function(x, ...) {
   if( "summary.tobit2" %in% class( x ) ) {
      model <- 2
   } else if( "summary.tobit5" %in% class( x ) ) {
      model <- 5
   } else {
      stop( "print.summary.selection currently works for tobit2 and tobit5",
         " models only" )
   }

   cat("--------------------------------------------\n")
   cat("Tobit ", model, " selection model/Maximum Likelihood estimation\n",
      sep = "" )
   cat(x$type, ", ", x$iterations, " iterations\n", sep="")
   cat("Return code ", x$code, ": ", x$message, "\n", sep="")
   if(!is.null(x$estimate)) {
      cat("Log-Likelihood:", loglikValue(x), "\n")
      if( model == 2 ) {
         cat(x$NObs, " observations (", x$N0, " censored and ", x$N1,
            " observed)", sep = "" )
      } else if( model == 5 ) {
         cat(x$NObs, " observations (", x$N1, " selection 1 and ", x$N2,
            " selection 2)", sep = "" )
      }
      cat( " and ", x$NActivePar, " free parameters (df = ",
          x$NObs - x$NActivePar, ")\n", sep="")
      cat("Probit selection equation:\n")
      printCoefmat(x$estimateS, signif.legend=FALSE)
      if( model == 2 ) {
         cat("Outcome equation:\n")
         printCoefmat(x$estimateO, signif.legend=FALSE)
      } else if( model == 5 ) {
         cat("Outcome equation 1:\n")
         printCoefmat(x$estimateO1, signif.legend=FALSE)
         cat("Outcome equation 2:\n")
         printCoefmat(x$estimateO2, signif.legend=FALSE)
      }
      cat("Error terms data:\n")
      printCoefmat(x$estimateErr)
   }
   cat("--------------------------------------------\n")
}
