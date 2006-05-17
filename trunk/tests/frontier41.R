library( micEcon )
data( Coelli )
Coelli$logOutput  <- log( Coelli$output )
Coelli$logCapital <- log( Coelli$capital )
Coelli$logLabour  <- log( Coelli$labour )

insFile <- file()
dtaFile  <- file()

writeFront41in( Coelli, "firm", "time", "logOutput",
   c( "logCapital", "logLabour" ), insFile = insFile, dtaFile = dtaFile  )

print( readLines( insFile ) )
print( readLines( dtaFile ) )
close( insFile )
close( dtaFile )




