heckit <- function( formula, probitformula, data ) {

   if( class( formula ) != "formula" ) {
      stop( "argument 'formula' must be a formula" )
   } else if( length( formula ) != 3 ) {
      stop( "argument 'formula' must be a 2-sided formula" )
   } else if( "probit" %in% substr( all.vars( formula ), 1, 6 ) ) {
      stop( paste( "argument 'formula' may not include variable names",
      "starting with 'probit'" ) )
   } else if( class( probitformula ) != "formula" ) {
      stop( "argument 'probitformula' must be a formula" )
   } else if( length( probitformula ) != 3 ) {
      stop( "argument 'probitformula' must be a 2-sided formula" )
   } else if( "probit" %in% substr( all.vars( probitformula ), 1, 6 ) ) {
      stop( paste( "argument 'probitformula' may not include a variable",
         "names starting with 'probit'" ) )
   }

   result <- list()

   data$probitdummy <- model.frame( probitformula, data = data )[ , 1 ]
   test <- levels( as.factor( as.numeric( data$probitdummy ) ) )
   if( length( test ) != 2 ) {
      stop( paste( "The left hand side of 'probitformula' may only contain",
         "1 and 0 or TRUE and FALSE" ) )
   } else if( !all.equal( test, c( "0", "1" ) ) ) {
      stop( paste( "The left hand side of 'probitformula' may only contain",
         "1 and 0 or TRUE and FALSE" ) )
   }

   result$probit <- glm( probitformula, binomial( link = "probit" ), data )

   data$probitlambda <- dnorm( result$probit$linear.predictors ) /
      pnorm( result$probit$linear.predictors )

   data$probitdelta <- data$probitlambda * ( data$probitlambda +
      result$probit$linear.predictors )

   step2formula <- as.formula( paste( formula[ 2 ], "~", formula[ 3 ],
      "+ probitlambda" ) )

   result$lm <- lm( step2formula, data, data$probitdummy == 1 )

   result$sigma <- sqrt( crossprod( residuals( result$lm ) ) /
      sum( data$probitdummy == 1 ) +
      mean( data$probitdelta[ data$probitdummy == 1 ] ) *
      coefficients( result$lm )[ "probitlambda" ]^2 )

   result$rho <- coefficients( result$lm )[ "probitlambda" ] / result$sigma

   class( result ) <- "heckit"
   return( result )
}

summary.heckit <- function( object, ... ) {
   cat( "1st step probit estimation:\n" )
   print( summary( object$probit, ... ) )
   cat( "\n2nd step linear estimation:\n" )
   print( summary( object$lm, ... ) )
   invisible( object )
}
