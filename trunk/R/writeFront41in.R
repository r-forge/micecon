writeFront41in <- function( data, crossSectionName, timePeriodName,
   yName, xNames = NULL, zNames = NULL,
   translog = FALSE, quadHalf = TRUE,
   functionType = 1, modelType = 1, logDepVar = "y", mu = "y", eta = "y",
   insFilename = "front41.ins", dtaFilename = sub( "\.ins$", ".dta", insFilename ),
   outFilename = sub( "\.ins$", ".out", insFilename ) ) {

   checkNames( c( crossSectionName, timePeriodName, yName, xNames, zNames ),
      names( data ) )

   if( !modelType %in% c( 1, 2 ) ) {
      stop( "Argument 'modelType' must be either 1 or 2." )
   }
   if( !functionType %in% c( 1, 2 ) ) {
      stop( "Argument 'functionType' must be either 1 or 2." )
   }
   if( !logDepVar %in% c( "y", "n" ) ) {
      stop( "Argument 'logDepVar' must be either 'y' or 'n'." )
   }
   if( !mu %in% c( "y", "n" ) ) {
      stop( "Argument 'mu' must be either 'y' or 'n'." )
   }
   if( modelType == 1 ) {
      if( !eta %in% c( "y", "n" ) ) {
         stop( "Argument 'eta' must be either 'y' or 'n'." )
      }
   }


   nCrossSection <- max( data[[ crossSectionName ]] )
   nTimePeriods  <- max( data[[ timePeriodName ]] )
   nTotalObs     <- nrow( data )
   nXvars        <- length( xNames )
   nXtotal       <- ifelse( translog, nXvars + nXvars * ( nXvars + 1 ) / 2,
                            nXvars )
   nZvars        <- length( zNames )

   if( modelType == 2 ) {
      eta <- nZvars
   }

   commentRow <- max( 16, nchar( dtaFilename ) + 1 )

   cat( modelType, rep( " ", commentRow - 1 ),
      "1=ERROR COMPONENTS MODEL, 2=TE EFFECTS MODEL\n",
      file = insFilename, sep = "" )
   cat( dtaFilename, rep( " ", commentRow - nchar( dtaFilename ) ),
      "DATA FILE NAME\n", file = insFilename, append = TRUE, sep = "" )
   cat( outFilename, rep( " ", commentRow - nchar( outFilename ) ),
      "OUTPUT FILE NAME\n", file = insFilename, append = TRUE, sep = "" )
   cat( functionType, rep( " ", commentRow - 1 ),
      "1=PRODUCTION FUNCTION, 2=COST FUNCTION\n",
      file = insFilename, append = TRUE, sep = "" )
   cat( logDepVar, rep( " ", commentRow - 1 ),
      "LOGGED DEPENDENT VARIABLE (Y/N)\n",
      file = insFilename, append = TRUE, sep = "" )
   cat( nCrossSection,
      rep( " ", commentRow - nchar( as.character( nCrossSection ) ) ),
      "NUMBER OF CROSS-SECTIONS\n",
      file = insFilename, append = TRUE, sep = "" )
   cat( nTimePeriods,
      rep( " ", commentRow - nchar( as.character( nTimePeriods ) ) ),
      "NUMBER OF TIME PERIODS\n",
      file = insFilename, append = TRUE, sep = "" )
   cat( nTotalObs,
      rep( " ", commentRow - nchar( as.character( nTotalObs ) ) ),
      "NUMBER OF OBSERVATIONS IN TOTAL\n",
      file = insFilename, append = TRUE, sep = "" )
   cat( nXtotal,
      rep( " ", commentRow - nchar( as.character( nXtotal ) ) ),
      "NUMBER OF REGRESSOR VARIABLES (Xs)\n",
      file = insFilename, append = TRUE, sep = "" )
   cat( mu, rep( " ", commentRow - 1 ),
      "MU (Y/N) [OR DELTA0 (Y/N) IF USING TE EFFECTS MODEL]\n",
      file = insFilename, append = TRUE, sep = "" )
   cat( eta, rep( " ", commentRow - nchar( as.character( eta ) ) ),
      "ETA (Y/N) [OR NUMBER OF TE EFFECTS REGRESSORS (Zs)]\n",
      file = insFilename, append = TRUE, sep = "" )
   cat( "n", rep( " ", commentRow - 1 ),
      "STARTING VALUES (Y/N)\n",
      file = insFilename, append = TRUE, sep = "" )

   dataTable <- cbind( data[[ crossSectionName ]], data[[ timePeriodName ]],
      data[[ yName ]] )

   if( nXvars > 0 ) {
      for( i in 1:nXvars ) {
         dataTable <- cbind( dataTable, data[[ xNames[ i ] ]] )
      }
      if( translog ) {
         for( i in 1:nXvars ) {
            for( j in i:nXvars ) {
               dataTable <- cbind( dataTable,
                  ifelse( i == j, 1 , 2 ) * ifelse( quadHalf, 0.5, 1 ) *
                  data[[ xNames[ i ] ]] * data[[ xNames[ j ] ]] )
            }
         }
      }
   }
   if( nZvars > 0 ) {
      for( i in 1:nZvars ) {
         dataTable <- cbind( dataTable, data[[ zNames[ i ] ]] )
      }
   }
   write.table( dataTable, file = dtaFilename, row.names = FALSE,
      col.names = FALSE, sep = "\t" )
}