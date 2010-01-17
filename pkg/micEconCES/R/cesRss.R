cesRss <- function( par, yName, xNames, data, vrs ) {

   yHat <- cesCalc( xNames = xNames, data = data, coef = par )

   result <- sum( ( data[[ yName ]] - yHat )^2 )
   if( is.na( result ) ) {
      result <- Inf
   }
   return( result )
}
