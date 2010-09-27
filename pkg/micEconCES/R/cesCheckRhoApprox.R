cesCheckRhoApprox <- function( rhoApprox, nElem ) {

   if( !is.vector( rhoApprox ) || length( rhoApprox ) != nElem ||
         !is.numeric( rhoApprox ) ) {
      stop( "argument 'rhoApprox' must be a numeric vector with exactly ",
         nElem, " elements" )
   }

   return( NULL )
}
