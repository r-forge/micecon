library( micEcon )

## *****************************
## Testing front41WriteInput

data( Coelli )
Coelli$logOutput  <- log( Coelli$output )
Coelli$logCapital <- log( Coelli$capital )
Coelli$logLabour  <- log( Coelli$labour )

insFile <- file()
dtaFile  <- file()

front41WriteInput( Coelli, "firm", "time", "logOutput",
   c( "logCapital", "logLabour" ), insFile = insFile, dtaFile = dtaFile  )

print( readLines( insFile ) )
print( readLines( dtaFile ) )

# irregular firm (cross section) identifier
set.seed( 20061705 )
Coelli$firm <- sample( c( 1:( nrow( Coelli ) + 20 ) ) )[ 1:nrow( Coelli ) ]

front41WriteInput( Coelli, "firm", "time", "logOutput",
   c( "logCapital", "logLabour" ), insFile = insFile, dtaFile = dtaFile  )

print( readLines( insFile ) )
print( readLines( dtaFile ) )

close( insFile )
close( dtaFile )


## *****************************
## Testing front41ReadOutput

outFile <- system.file( "front41/EG1.OUT", package = "micEcon" )
sfa <- front41ReadOutput( outFile )
print( sfa )
