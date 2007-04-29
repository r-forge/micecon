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
estResultLA <- summary( aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86[ set, ], method = "LA:SL",
   maxiter = 1, rcovformula = 1, tol = 1e-7 ), elaFormula = "Ch", quantNames = wNames )
print( estResultLA )
# imposing restrictions via TX
estResultLATX <- summary( aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86[ set, ], method = "LA:SL",
   maxiter = 1, rcovformula = 1, tol = 1e-7, TX = TRUE ), elaFormula = "Ch", quantNames = wNames )
print( estResultLATX )
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
estResultLAhom <- summary(  aidsEst( pNames, wNames, "xFood", sym = FALSE,
   data = Blanciforti86[ set, ], method = "LA:SL",
   maxiter = 1, rcovformula = 1, tol = 1e-7 ), elaFormula = "Ch", quantNames = wNames )
print( estResultLAhom )
# imposing restrictions via TX
estResultLAhomTX <- summary(  aidsEst( pNames, wNames, "xFood", sym = FALSE,
   data = Blanciforti86[ set, ], method = "LA:SL",
   maxiter = 1, rcovformula = 1, tol = 1e-7, TX = TRUE ), elaFormula = "Ch", quantNames = wNames )
print( estResultLAhomTX )
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
estResultLAunr <- summary( aidsEst( pNames, wNames, "xFood", hom = FALSE, sym = FALSE,
   data = Blanciforti86[ set, ], method = "LA:SL",
   maxiter = 1, rcovformula = 1, tol = 1e-7 ), elaFormula = "Ch", quantNames = wNames )
print( estResultLAunr )
# imposing restrictions via TX
estResultLAunrTX <- summary( aidsEst( pNames, wNames, "xFood", hom = FALSE, sym = FALSE,
   data = Blanciforti86[ set, ], method = "LA:SL",
   maxiter = 1, rcovformula = 1, tol = 1e-7, TX = TRUE ), elaFormula = "Ch", quantNames = wNames )
print( estResultLAunrTX )
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
estResultAIDS <- summary( aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86[ setWo1, ], maxiter = 1,
   rcovformula=1, tol=1e-7,
   method = "IL:L" ), elaFormula = "AIDS", quantNames = wNames )
print( estResultAIDS )
# imposing restrictions via TX
estResultAIDSTX <- summary( aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86[ setWo1, ], maxiter = 1,
   rcovformula=1, tol=1e-7,
   method = "IL:L", TX = TRUE ), elaFormula = "AIDS", quantNames = wNames )
print( estResultAIDSTX )
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
estResultAIDShom <- summary( aidsEst( pNames, wNames, "xFood", sym = FALSE,
   data = Blanciforti86[ setWo1, ], maxiter = 1,
   rcovformula=1, tol=1e-7,
   method = "IL:L" ), elaFormula = "AIDS", quantNames = wNames )
print( estResultAIDShom )
# imposing restrictions via TX
estResultAIDShomTX <- summary( aidsEst( pNames, wNames, "xFood", sym = FALSE,
   data = Blanciforti86[ setWo1, ], maxiter = 1,
   rcovformula=1, tol=1e-7,
   method = "IL:L", TX = TRUE ), elaFormula = "AIDS", quantNames = wNames )
print( estResultAIDShomTX )
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
estResultAIDSunr <- summary( aidsEst( pNames, wNames, "xFood", hom = FALSE, sym = FALSE,
   data = Blanciforti86[ setWo1, ], maxiter = 1,
   rcovformula=1, tol=1e-7,
   method = "IL:L" ), elaFormula = "AIDS", quantNames = wNames )
print( estResultAIDSunr )
# imposing restrictions via TX
estResultAIDSunrTX <- summary( aidsEst( pNames, wNames, "xFood", hom = FALSE, sym = FALSE,
   data = Blanciforti86[ setWo1, ], maxiter = 1,
   rcovformula=1, tol=1e-7,
   method = "IL:L", TX = TRUE ), elaFormula = "AIDS", quantNames = wNames )
print( estResultAIDSunrTX )
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
estResultLaSNa <- summary( aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86, maxiter = 1,
   rcovformula=1, tol=1e-7, method = "LA:S" ), elaFormula = "AIDS", quantNames = wNames )
print( estResultLaSNa )

estResultLaSlNa <- summary( aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86, maxiter = 1,
   rcovformula=1, tol=1e-7, method = "LA:SL" ), elaFormula = "AIDS", quantNames = wNames )
print( estResultLaSlNa )

estResultLaLNa <- summary( aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86, maxiter = 1,
   rcovformula=1, tol=1e-7, method = "LA:L" ), elaFormula = "AIDS", quantNames = wNames )
print( estResultLaLNa )

estResultAIDSNa <- summary( aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86, maxiter = 1,
   rcovformula=1, tol=1e-7, method = "IL:L" ), elaFormula = "AIDS", quantNames = wNames )
print( estResultAIDSNa )


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

