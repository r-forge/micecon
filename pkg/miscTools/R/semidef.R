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
      } else if( !isSymmetric( m ) ) {
         stop( "argument 'm' must be a symmetric matrix" )
      }
      n <- nrow( m )
      if( !positive ) {
         m <- -m
      }
      if( method == "det" ) {
         result <- TRUE
         for( i in 1:n ) {
            comb <- combn( n, i )
            for( j in 1:ncol( comb ) ) {
               princMin <-  det( m[ comb[ , j ], comb[ , j ], drop = FALSE ] )
               result <- result && ( princMin >= -tol )
            }
         }
      } else if( method == "eigen" ) {
         ev <- eigen( m, only.values = TRUE )$values
         if( is.complex( ev ) ) {
            stop( "complex (non-real) eigenvalues,",
               " which could be caused by a non-symmetric matrix" )
         }
         result <- min( ev ) > -tol
      } else {
         stop( "argument 'method' must be either 'det' or 'eigen'" )
      }
   }
   return( result )
}
