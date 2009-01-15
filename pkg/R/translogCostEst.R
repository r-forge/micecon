translogCostEst <- function( cName, yName, pNames, data, 
   fNames = NULL, shifterNames = NULL,
   dataLogged = FALSE, homPrice = TRUE, ... ) {

   checkNames( c( cName, yName, pNames, fNames, shifterNames ), names( data ) )

   if( dataLogged ) {
      logData   <- data
   } else {
      if( "plm.dim" %in% class( data ) ) {
         logData <- data[ , 1:2 ]
      } else {
         logData <- data.frame( no = c( 1:nrow( data ) ) )
      }
      logData[[ cName ]] <- log( data[[ cName ]] )
      logData[[ yName ]] <- log( data[[ yName ]] )
      for( i in seq( along = pNames ) ) {
         logData[[ pNames[ i ] ]] <- log( data[[ pNames[ i ] ]] )
      }
      for( i in seq( along = fNames ) ) {
         logData[[ fNames[ i ] ]] <-
            log( data[[ fNames[ i ] ]] )
      }
      for( i in seq( along = shifterNames ) ) {
         logData[[ shifterNames[ i ] ]] <-
            log( data[[ shifterNames[ i ] ]] )
      }
   }

   qxVars <- yName
   if( homPrice ){
      stop( "imposing linear homogeneity in input prices has not been",
         " implemented yet" )
   } else {
      qxVars <- c( qxVars, pNames )
   }
   qxVars <- c( qxVars, fNames )
      
   result <- quadFuncEst( yName = cName, xNames = qxVars, 
      data = logData, shifterNames = shifterNames, quadHalf = TRUE, ... )

   result$r2nonLog <- rSquared( exp( logData[[ yName ]] ),
      exp( logData[[ yName ]] ) - exp( result$fitted ) )

   if( !dataLogged ){
      result$fitted <- exp( result$fitted )
   }

   class( result ) <- "translogCostEst"
   return( result )
}
