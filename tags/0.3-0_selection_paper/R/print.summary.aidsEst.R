print.summary.aidsEst <- function( x, ... ) {
   print.aidsEst( x )
   cat( "\nDemand Elasticities " )
   if( x$ela$formula == "Ch" ) {
      cat( "(Formula of Chalfant / Goddard)\n" )
   } else if( x$ela$formula == "AIDS" ) {
      cat( "(original AIDS formula)\n" )
   }
   cat( "Expenditure Elasticities\n" )
   print( x$ela$exp )
   cat( "\nMarshallian (uncompensated) Price Elasticities\n" )
   print( x$ela$marshall )
   cat( "\nHicksian (compensated) Price Elasticities\n" )
   print( x$ela$hicks )
   invisible( x )
}
