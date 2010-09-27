cesCheckRhoApprox <- function( rhoApprox, elemNames ) {

   if( !is.vector( rhoApprox ) || !is.numeric( rhoApprox ) ) {
      stop( "argument 'rhoApprox' must be a numeric vector" )
   } else if( is.null( names( rhoApprox ) ) ) {
      stop( "elements of argument 'rhoApprox' must be named" )
   }
   elemInRhoApprox <- elemNames %in% names( rhoApprox )
   if( !all( elemInRhoApprox ) ) {
      stop( "following elements are missing in argument 'rhoApprox': ",
         paste( elemNames[ !elemInRhoApprox ], collapse = ", " ) )
   }

   return( NULL )
}
