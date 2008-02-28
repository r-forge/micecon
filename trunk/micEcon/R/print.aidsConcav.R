print.aidsConcav <- function( x, ... ) {
   cat( "\nChecking concavity of an estimated " )
   cat( "Almost Ideal Demand System (AIDS):\n" )
   cat( "Concavity is fulfilled at " )
   cat( x$nConcavObs, "out of", x$nValidObs, "observations" )
   cat( " (", x$concavPercent, "%)\n", sep = "" )
}
