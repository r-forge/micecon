library( micEcon )
data( Blanciforti86 )
options( digits = 3 )

set <- !is.na( Blanciforti86$pFood1 )
setWo1 <- set & rownames( Blanciforti86 ) != 1947
pNames <- c( "pFood1", "pFood2", "pFood3", "pFood4" )
wNames <- c( "wFood1", "wFood2", "wFood3", "wFood4" )
wMeans <- colMeans( Blanciforti86[ set, wNames ] )
pMeans <- colMeans( Blanciforti86[ set, pNames ] )

cat( paste( "\nRepeating the demand analysis of Blanciforti, Green",
   "& King (1986)\n" ) )
estResultLA <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86[ set, ], method = "LA:SL" )
print( estResultLA )
print( summary( estResultLA, elaFormula = "Ch", quantNames = wNames ) )
# imposing restrictions via TX
estResultLATX <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86[ set, ], method = "LA:SL",
   TX = TRUE )
print( estResultLATX )
print( summary( estResultLATX, elaFormula = "Ch", quantNames = wNames ) )
estResultLATX$call <- NULL
estResultLATX$est$bt <- NULL
estResultLATX$est$btcov <- NULL
estResultLATX$est$x <- NULL
estResultLATX$est$TX <- NULL
estResultLATX$est$q.restr <- NULL
estResultLA$call <- NULL
estResultLA$est$x <- NULL
estResultLA$est$R.restr <- NULL
estResultLA$est$q.restr <- NULL
print( all.equal( estResultLA, estResultLATX ) )

## only homogeneity (no symmetry imposed)
estResultLAhom <- aidsEst( pNames, wNames, "xFood", sym = FALSE,
   data = Blanciforti86[ set, ], method = "LA:SL" )
print( estResultLAhom )
print( summary( estResultLAhom, elaFormula = "Ch", quantNames = wNames ) )
# imposing restrictions via TX
estResultLAhomTX <- aidsEst( pNames, wNames, "xFood", sym = FALSE,
   data = Blanciforti86[ set, ], method = "LA:SL",
   TX = TRUE )
print( estResultLAhomTX )
print( summary( estResultLAhomTX, elaFormula = "Ch", quantNames = wNames ) )
estResultLAhomTX$call <- NULL
estResultLAhomTX$est$bt <- NULL
estResultLAhomTX$est$btcov <- NULL
estResultLAhomTX$est$x <- NULL
estResultLAhomTX$est$TX <- NULL
estResultLAhomTX$est$q.restr <- NULL
estResultLAhom$call <- NULL
estResultLAhom$est$x <- NULL
estResultLAhom$est$R.restr <- NULL
estResultLAhom$est$q.restr <- NULL
print( all.equal( estResultLAhom, estResultLAhomTX ) )

## unrestricted (no homogeneity and no symmetry imposed)
estResultLAunr <- aidsEst( pNames, wNames, "xFood", hom = FALSE, sym = FALSE,
   data = Blanciforti86[ set, ], method = "LA:SL" )
print( estResultLAunr )
print( summary( estResultLAunr, elaFormula = "Ch", quantNames = wNames ) )
# imposing restrictions via TX
estResultLAunrTX <- aidsEst( pNames, wNames, "xFood", hom = FALSE, sym = FALSE,
   data = Blanciforti86[ set, ], method = "LA:SL",
   TX = TRUE )
print( estResultLAunrTX )
print( summary( estResultLAunrTX, elaFormula = "Ch", quantNames = wNames ) )
estResultLAunrTX$call <- NULL
estResultLAunrTX$est$bt <- NULL
estResultLAunrTX$est$btcov <- NULL
estResultLAunrTX$est$x <- NULL
estResultLAunrTX$est$TX <- NULL
estResultLAunrTX$est$q.restr <- NULL
estResultLAunr$call <- NULL
estResultLAunr$est$x <- NULL
estResultLAunr$est$R.restr <- NULL
estResultLAunr$est$q.restr <- NULL
print( all.equal( estResultLAunr, estResultLAunrTX ) )

#####################################################
cat( paste( "\nRepeating the evaluation of different elasticity formulas",
   "of Green & Alston (1990): iterated AIDS\n" ) )
estResultAIDS <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86[ setWo1, ], method = "IL:L" )
print( estResultAIDS )
print( summary( estResultAIDS, elaFormula = "AIDS", quantNames = wNames ) )
# imposing restrictions via TX
estResultAIDSTX <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86[ setWo1, ], method = "IL:L", TX = TRUE )
print( estResultAIDSTX )
print( summary( estResultAIDSTX, elaFormula = "AIDS", quantNames = wNames ) )
estResultAIDSTX$call <- NULL
estResultAIDSTX$est$bt <- NULL
estResultAIDSTX$est$btcov <- NULL
estResultAIDSTX$est$x <- NULL
estResultAIDSTX$est$TX <- NULL
estResultAIDSTX$est$q.restr <- NULL
estResultAIDS$call <- NULL
estResultAIDS$est$x <- NULL
estResultAIDS$est$R.restr <- NULL
estResultAIDS$est$q.restr <- NULL
print( all.equal( estResultAIDS, estResultAIDSTX ) )

## only homogeneity (no symmetry imposed)
estResultAIDShom <- aidsEst( pNames, wNames, "xFood", sym = FALSE,
   data = Blanciforti86[ setWo1, ], method = "IL:L" )
print( estResultAIDShom )
print( summary( estResultAIDShom, elaFormula = "AIDS", quantNames = wNames ) )
# imposing restrictions via TX
estResultAIDShomTX <- aidsEst( pNames, wNames, "xFood", sym = FALSE,
   data = Blanciforti86[ setWo1, ], method = "IL:L", TX = TRUE )
print( estResultAIDShomTX )
print( summary( estResultAIDShomTX, elaFormula = "AIDS", quantNames = wNames ) )
estResultAIDShomTX$call <- NULL
estResultAIDShomTX$est$bt <- NULL
estResultAIDShomTX$est$btcov <- NULL
estResultAIDShomTX$est$x <- NULL
estResultAIDShomTX$est$TX <- NULL
estResultAIDShomTX$est$q.restr <- NULL
estResultAIDShom$call <- NULL
estResultAIDShom$est$x <- NULL
estResultAIDShom$est$R.restr <- NULL
estResultAIDShom$est$q.restr <- NULL
print( all.equal( estResultAIDShom, estResultAIDShomTX ) )

## unrestricted (no homogeneity and no symmetry imposed)
estResultAIDSunr <- aidsEst( pNames, wNames, "xFood", hom = FALSE, sym = FALSE,
   data = Blanciforti86[ setWo1, ], method = "IL:L" )
print( estResultAIDSunr )
print( summary( estResultAIDSunr, elaFormula = "AIDS", quantNames = wNames ) )
# imposing restrictions via TX
estResultAIDSunrTX <- aidsEst( pNames, wNames, "xFood", hom = FALSE, sym = FALSE,
   data = Blanciforti86[ setWo1, ], method = "IL:L", TX = TRUE )
print( estResultAIDSunrTX )
print( summary( estResultAIDSunrTX, elaFormula = "AIDS", quantNames = wNames ) )
estResultAIDSunrTX$call <- NULL
estResultAIDSunrTX$est$bt <- NULL
estResultAIDSunrTX$est$btcov <- NULL
estResultAIDSunrTX$est$x <- NULL
estResultAIDSunrTX$est$TX <- NULL
estResultAIDSunrTX$est$q.restr <- NULL
estResultAIDSunr$call <- NULL
estResultAIDSunr$est$x <- NULL
estResultAIDSunr$est$R.restr <- NULL
estResultAIDSunr$est$q.restr <- NULL
print( all.equal( estResultAIDSunr, estResultAIDSunrTX ) )

## with NAs
estResultLaSNa <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86,
   method = "LA:S" )
print( estResultLaSNa )
print( summary( estResultLaSNa, elaFormula = "AIDS", quantNames = wNames ) )

estResultLaSlNa <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86,
   method = "LA:SL" )
print( estResultLaSlNa )
print( summary( estResultLaSlNa, elaFormula = "AIDS", quantNames = wNames ) )

estResultLaLNa <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86,
   method = "LA:L" )
print( estResultLaLNa )
print( summary( estResultLaLNa, elaFormula = "AIDS", quantNames = wNames ) )

estResultAIDSNa <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86, method = "IL:L" )
print( estResultAIDSNa )
print( summary( estResultAIDSNa, elaFormula = "AIDS", quantNames = wNames ) )


########## Elasticities ###############
cat( "\nAIDS: Elasticities\n" )
wMeans <- colMeans( Blanciforti86[ set, wNames ] )
ela <- aidsEla( estResultAIDS$coef, wMeans, pMeans, formula = "AIDS" )
print( ela )
elaTX <- aidsEla( estResultAIDSTX$coef, wMeans, pMeans, formula = "AIDS" )
print( elaTX )
print( all.equal( ela, elaTX ) )


cat( "\nLA: Elasticity formula of non-linear AIDS\n" )
wMeans <- colMeans( Blanciforti86[ set, wNames ] )
ela <- aidsEla( estResultLA$coef, wMeans, pMeans, formula = "AIDS" )
print( ela )
elaTX <- aidsEla( estResultLATX$coef, wMeans, pMeans, formula = "AIDS" )
print( elaTX )
print( all.equal( ela, elaTX ) )

cat( "\n********** Elasticities ***************" )
cat( "\nLA: Elasticity formula of Chalfant / Goddard\n" )
ela <- aidsEla( estResultLA$coef, wMeans, formula = "Ch" )
print( ela )

cat( "\nLA: Elasticity formula of Eales + Unnevehr\n" )
wMeans <- colMeans( Blanciforti86[ set, wNames ] )
ela <- aidsEla( estResultLA$coef, wMeans, formula = "EU" )
print( ela )


############# Price indices ##############
options( digits = 5 )
cat( "\n************** Price indices **************\n" )
cat( "\nStone index\n" )
pxS <- aidsPx( "S", pNames, wNames, Blanciforti86 )
print( pxS )

cat( "\nStone index with lagged shares\n" )
pxSL <- aidsPx( "SL", pNames, wNames, Blanciforti86 )
print( pxSL )

cat( "\nPaasche index\n" )
pxP <- aidsPx( "P", pNames, wNames, Blanciforti86 )
print( pxP )

cat( "\nLaspeyres index\n" )
pxL <- aidsPx( "L", pNames, wNames, Blanciforti86 )
print( pxL )

cat( "\nTornqvist index\n" )
pxT <- aidsPx( "T", pNames, wNames, Blanciforti86 )
print( pxT )

cat( "\nTranslog index\n" )
pxTL <- aidsPx( "TL", pNames, wNames, Blanciforti86,
   coef = estResultLA$coef )
print( pxTL )

########### fitted values #################
options( digits = 3 )
fittedAIDS <- aidsCalc( pNames, "xFood", Blanciforti86[ -1, ],
   coef = estResultAIDS$coef )
print( fittedAIDS )
if( max( abs( fittedAIDS$shares[ !is.na( fittedAIDS$shares ) ] -
   estResultAIDS$wFitted ) ) > 1e-5 ) {
   stop( "fitted shares of AIDS are wrong" )
}
if( max( abs( fittedAIDS$quant[ !is.na( fittedAIDS$quant ) ] -
   estResultAIDS$qFitted ) ) > 1e-5 ) {
   stop( "fitted quantities of AIDS are wrong" )
}
fittedAIDSTX <- aidsCalc( pNames, "xFood", Blanciforti86[ -1, ],
   coef = estResultAIDSTX$coef )
print( fittedAIDSTX )
print( all.equal( fittedAIDS, fittedAIDSTX ) )

fittedLA <- aidsCalc( pNames, "xFood", Blanciforti86[ set, ],
   coef = estResultLA$coef, lnp = estResultLA$lnp )
print( fittedLA )
if( max( abs( fittedLA$shares[ -1, ] - estResultLA$wFitted[ setWo1, ] ),
      na.rm = TRUE ) > 1e-5 ) {
   stop( "fitted shares of LA-AIDS are wrong" )
}
if( max( abs( fittedLA$quant[ -1, ] - estResultLA$qFitted[ setWo1, ] ),
      na.rm = TRUE ) > 1e-5 ) {
   stop( "fitted quantities of LA-AIDS are wrong" )
}
fittedLATX <- aidsCalc( pNames, "xFood", Blanciforti86[ set, ],
   coef = estResultLATX$coef, lnp = estResultLATX$lnp )
print( fittedLATX )
print( all.equal( fittedLA, fittedLATX ) )

####### consistency ###################
consist <- aidsTestConsist( pNames, wNames, "xFood", Blanciforti86[ set, ],
   coef = estResultAIDS$coef )
print( consist )
class( consist ) <- NULL
print( consist )

