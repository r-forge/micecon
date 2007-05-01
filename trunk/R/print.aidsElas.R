print.aidsElas <- function( x, ... ) {

   cat( "\nDemand Elasticities " )
   if( x$method == "Ch" ) {
      cat( "(formula of Chalfant / Goddard)\n" )
   } else if( x$method == "EU" ) {
      cat( "(formula of Eales and Unnevehr)\n" )
   } else if( x$method == "AIDS" ) {
      cat( "(original AIDS formula)\n" )
   } else {
      cat( "(unknown formula '", x$method, "')\n", sep = "" )
   }
   cat( "Expenditure Elasticities\n" )
   print( x$exp )
   cat( "\nMarshallian (uncompensated) Price Elasticities\n" )
   print( x$marshall )
   cat( "\nHicksian (compensated) Price Elasticities\n" )
   print( x$hicks )
   invisible( x )
}
