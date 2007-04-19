summary.selection <- function(object, ...) {
   ## wrapper function for "summary.tobit" and "summary.heckit"
   ## object      object of class "selection"
   ## ...         additional arguments for "summary.tobit" and "summary.heckit"

   if( object$method == "ml" ) {
      result <- summary.tobit( object, ... )
   } else if( object$method == "2step" )  {
      result <- summary.heckit( object, ... )
   }
   result$method <- object$method
   class( result ) <- c( "summary.selection", class( result ) )
   return( result )
}

print.summary.selection <- function(x, ...) {

   cat("--------------------------------------------\n")
   if( x$method == "ml" ) {
      print.summary.tobit( x, ... )
   } else if( x$method == "2step" ) {
      print.summary.heckit( x, ... )
   }
   cat("--------------------------------------------\n")
   invisible( x )
}
