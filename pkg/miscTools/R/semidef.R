isSemidefinite <- function( m, ... )
   UseMethod( "isSemidefinite" )

isSemidefinite.default <- function( m, ... ) {
   stop( "there is currently no default method available" )
}

## ----- test positive / negative semidefiniteness
isSemidefinite.matrix <- function( m, positive = TRUE,
   tol = 100 * .Machine$double.eps,
   method = ifelse( nrow( m ) < 13, "det", "eigen" ), ... ) {

   if( !is.matrix( m ) ) {
      stop( "argument 'm' must be a matrix" )
   } else {
      if( nrow( m ) != ncol( m ) ) {
         stop( "argument 'm' or each of its elements must be a _quadratic_ matrix" )
      } else if( !isSymmetric( unname( m ), tol = 1000 * tol ) ) {
         stop( "argument 'm' must be a symmetric matrix" )
      }
      # make sure that the matrix is almost exactly symmetric
      # even if it is slightly non-symmetric 
      m <- ( m + t(m) ) / 2
      n <- nrow( m )
      if( !positive ) {
         m <- -m
      }
      if( n >= 12 && method == "det" ) {
         warning( "using method 'det' could take a very long time",
            " for matrices with more than 12 rows and columns;",
            " it is suggested to use method 'eigen' for larger matrices",
            immediate. = TRUE )
      }
      if( method == "det" ) {
         result <- TRUE
         for( i in 1:n ) {
            comb <- combn( n, i )
            for( j in 1:ncol( comb ) ) {
               mat <- m[ comb[ , j ], comb[ , j ], drop = FALSE ]
               if( rcond( mat ) >= tol ) {
                  princMin <-  det( mat )
               } else {
                  princMin <- 0
               }
               result <- result && ( princMin >= -tol )
            }
         }
      } else if( method == "eigen" ) {
         if( rcond( m ) >= tol || n == 1 ) {
            ev <- eigen( m, only.values = TRUE )$values
            if( is.complex( ev ) ) {
               stop( "complex (non-real) eigenvalues,",
                  " which could be caused by a non-symmetric matrix" )
            }
            result <- min( ev ) > -tol
         } else {
            result <- TRUE
            for( i in 1:n ) {
               mm <- m[ -i, -i, drop = FALSE ]
               result <- result &
                  semidefiniteness( mm, tol = tol, method = method  )
            }
         }
      } else {
         stop( "argument 'method' must be either 'det' or 'eigen'" )
      }
   }
   return( result )
}

isSemidefinite.list <- function( m, ... ) {
   
   if( !is.list( m ) ) {
      stop( "argument 'm' must be a matrix or a list of matrices" )
   } else if( !all( sapply( m, is.matrix ) ) ) {
      stop( "all components of the list specified by argument 'm'",
         " must be matrices" )
   }
      
   result <- logical( length( m ) )
   for( t in 1:length( m ) ) {
      result[ t ] <- isSemidefinite( m[[ t ]], ... )
   }
   return( result )
}

semidefiniteness <- function( m, ... ) {
   result <- isSemidefinite( m = m, ... )
   return( result )
}
