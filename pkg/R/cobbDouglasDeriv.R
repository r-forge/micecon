cobbDouglasDeriv <- function( xNames, data, coef,
   yName = NULL, dataLogged = FALSE ) {

   checkNames( c( xNames, yName ), names( data ) )

   nExog <- length( xNames )

   if( nExog + 1 != length( coef ) ) {
      stop( "a Cobb-Douglas function with ", nExog, " exogenous variables",
         " must have exactly ", nExog + 1, " coefficients" )
   }

   coefNames <- paste( "a", c( 0:nExog ), sep = "_" )
   if( is.null( names( coef ) ) ) {
      names( coef ) <- coefNames
   } else {
      coefMissing <- !( coefNames %in% names( coef ) )
      if( any( coefMissing ) ) {
         stop( "coefficient(s) ",
            paste( coefNames[ coefMissing ], collapse = ", " ),
            " are missing" )
      }
      rm( coefMissing )
   }
   rm( coefNames )

   if( dataLogged ) {
      logData <- data
   } else {
      logData <- .micEconLogData( data = data, varNames = xNames )
   }

   result <- list()

   if( is.null( yName ) ){
      logyHat <- cobbDouglasCalc( xNames = xNames, data = logData,
         coef = coef, dataLogged = TRUE )
   } else {
      if( dataLogged ) {
         logyHat <- data[[ yName ]]
      } else {
         logyHat <- log( data[[ yName ]] )
      }
   }

   deriv <- matrix( NA, nrow( data ), nExog )
   for( i in seq( along = xNames ) ) {
      deriv[ , i ] <- coef[ paste( "a", i, sep = "_" ) ] *
         exp( logyHat ) / exp( logData[[ xNames[ i ] ]] )
   }

   colnames( deriv ) <- xNames
   result$deriv      <- as.data.frame( deriv )

   class( result ) <- "cobbDouglasDeriv"
   return( result )
}
