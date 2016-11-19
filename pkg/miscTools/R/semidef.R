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
      if( method == "det" ) {
         if( positive ) {
            result <- ( min( diag( m ) ) >= -tol )
         } else {
            result <- ( max( diag( m ) ) <= tol )
         }
         if( n > 1 ) {
            for( i in 2:n ) {
               if( positive ) {
                  result <- result && ( det( m[ 1:i, 1:i ] ) >= -tol )
               } else {
                  result <- result && ( det( m[ 1:i, 1:i ] ) * ( -1 )^i >= -tol )
               }
            }
         }
      } else if( method == "eigen" ) {
         if( positive ) {
            result <- ( min( eigen( m )$values ) > -tol )
         } else {
            result <- ( max( eigen( m )$values ) < tol )
         }
      } else {
         stop( "argument 'method' must be either 'det' or 'eigen'" )
      }
   }
   return( result )
}
