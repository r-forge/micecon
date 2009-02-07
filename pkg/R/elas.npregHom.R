elas.npregHom <- function( object, data = NULL, yObs = FALSE, ... ) {

   if( is.null( data ) ) {
      data <- eval( object$call$data )
   }

   if( yObs ) {
      yValues = data[[ object$yName ]]
   } else {
      yValues <- fitted( object$est )
   }

   result <- as.data.frame( matrix( NA, nrow = nrow( object$grad ),
      ncol = ncol( object$grad ) ) )
   names( result ) <- colnames( object$grad )

   for( i in colnames( object$grad ) ) {
      result[[ i ]] <- object$grad[ , i ] * data[[ i ]] / yValues
   }

   return( result )
}