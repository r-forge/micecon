cesRss <- function( par, data ) {
   yHat <- cesCalc( xNames = names( data )[ -1 ], data = data, coef = par )
   return( sum( ( data$y - yHat )^2 ) )
}
