path <- "micEcon/R/"
micEconFiles <- list.files( path, all.files = TRUE )
micEconFiles <- grep( ".*\\.R$", micEconFiles, value = TRUE )
for( i in seq( along = micEconFiles ) ) {
   returnedByTry <- try( source( paste( path, micEconFiles[ i ], sep = "" ) ) )
   if( class( returnedByTry ) == "try-error" ) {
      cat( paste( micEconFiles[ i ], ": ", returnedByTry, sep = "" ) )
   }
}
