## ----- test positive / negative semidefiniteness
semidefiniteness <- function( m, positive = TRUE, tol = .Machine$double.eps,
      method = "det" ) {

   result <- list ()
   if( is.list( m ) ) {
      result <- logical( length( m ) )
      for( t in 1:length( m ) ) {
         result[ t ] <- semidefiniteness( m[[ t ]], positive = positive,
            tol = tol, method = method )
      }
   } else if( !is.matrix( m ) ) {
      stop( "argument 'm' must be a matrix or a list of matrices" )
   } else {
      if( nrow( m ) != ncol( m ) ) {
         stop( "argument 'm' or each of its elements must be a _quadratic_ matrix" )
      }
      n <- nrow( m )
      if( !positive ) {
         m <- -m
      }
      if( method == "det" ) {
         result <- ( min( diag( m ) ) >= -tol )
         if( n > 1 ) {
            for( i in 2:n ) {
               result <- result && ( det( m[ 1:i, 1:i ] ) >= -tol )
            }
         }
      } else if( method == "eigen" ) {
         result <- ( min( eigen( m, only.values = TRUE )$values ) > -tol )
      } else {
         stop( "argument 'method' must be either 'det' or 'eigen'" )
      }
   }
   return( result )
}
