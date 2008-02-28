print.aidsMono <- function( x, ... ) {
   cat( "\nChecking monotonicity of an estimated " )
   cat( "Almost Ideal Demand System (AIDS):\n" )
   cat( "Monotonicity is fulfilled at " )
   cat( x$nMonoObs, "out of", x$nValidObs, "observations" )
   cat( " (", x$monoPercent, "%)\n", sep = "" )
}
