cesRss <- function( par, data ) {

   yHat <- cesCalc( xNames = names( data )[ -1 ], data = data, coef = par )

   result <- sum( ( data$y - yHat )^2 )
   if( is.na( result ) ) {
      result <- Inf
   }
   return( result )
}
