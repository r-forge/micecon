library( micEcon )
data( Blanciforti86 )
options( digits = 3 )

set <- !is.na( Blanciforti86$pFood1 )
setWo1 <- set & rownames( Blanciforti86 ) != 1947
pNames <- c( "pFood1", "pFood2", "pFood3", "pFood4" )
wNames <- c( "wFood1", "wFood2", "wFood3", "wFood4" )
wMeans <- colMeans( Blanciforti86[ set, wNames ] )
pMeans <- colMeans( Blanciforti86[ set, pNames ] )


## estimations with different price indices
# AIDS: translog
estResultTl <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86[ set, ], method = "IL" )
print( estResultTl )
print( summary( estResultTl ) )

# LA-AIDS: Stone
estResultLaS <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86[ set, ], priceIndex = "S" )
print( estResultLaS )
print( summary( estResultLaS ) )

# LA-AIDS: Stone with lagged shares
estResultLaSl <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86[ set, ], priceIndex = "SL" )
print( estResultLaSl )
print( summary( estResultLaSl ) )

# LA-AIDS: Paasche
estResultLaP <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86[ set, ], priceIndex = "P" )
print( estResultLaP )
print( summary( estResultLaP ) )

# LA-AIDS: Laspeyres, simplified
estResultLaLs <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86[ set, ], priceIndex = "Ls" )
print( estResultLaLs )
print( summary( estResultLaLs ) )

# LA-AIDS: Tornqvist
estResultLaT <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86[ set, ], priceIndex = "T" )
print( estResultLaT )
print( summary( estResultLaT ) )


cat( paste( "\nRepeating the demand analysis of Blanciforti, Green",
   "& King (1986)\n" ) )
estResultLA <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86[ set, ], priceIndex = "SL" )
print( estResultLA )
print( summary( estResultLA ) )
print( elas( estResultLA, method = "Ch", quantNames = wNames ) )
# imposing restrictions via restrict.regMat
estResultLATX <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86[ set, ], priceIndex = "SL",
   restrict.regMat = TRUE )
print( estResultLATX )
print( summary( estResultLATX ) )
print( elas( estResultLATX, method = "Ch", quantNames = wNames ) )

## only homogeneity (no symmetry imposed)
estResultLAhom <- aidsEst( pNames, wNames, "xFood", sym = FALSE,
   data = Blanciforti86[ set, ], priceIndex = "SL" )
print( estResultLAhom )
print( summary( estResultLAhom ) )
print( elas( estResultLAhom, method = "Ch", quantNames = wNames ) )
# imposing restrictions via restrict.regMat
estResultLAhomTX <- aidsEst( pNames, wNames, "xFood", sym = FALSE,
   data = Blanciforti86[ set, ], priceIndex = "SL",
   restrict.regMat = TRUE )
print( estResultLAhomTX )
print( summary( estResultLAhomTX ) )
print( elas( estResultLAhomTX, method = "Ch", quantNames = wNames ) )

## unrestricted (no homogeneity and no symmetry imposed)
estResultLAunr <- aidsEst( pNames, wNames, "xFood", hom = FALSE, sym = FALSE,
   data = Blanciforti86[ set, ], priceIndex = "SL" )
print( estResultLAunr )
print( summary( estResultLAunr ) )
print( elas( estResultLAunr, method = "Ch", quantNames = wNames ) )
# imposing restrictions via restrict.regMat
estResultLAunrTX <- aidsEst( pNames, wNames, "xFood", hom = FALSE, sym = FALSE,
   data = Blanciforti86[ set, ], priceIndex = "SL",
   restrict.regMat = TRUE )
print( estResultLAunrTX )
print( summary( estResultLAunrTX ) )
print( elas( estResultLAunrTX, method = "Ch", quantNames = wNames ) )


#####################################################
## Estimation with demand shifters
Blanciforti86$trend <- c( 0:( nrow( Blanciforti86 ) - 1 ) )
estResultLAtrend <- aidsEst( pNames, wNames, "xFood",
   shifterNames = c( "trend" ), data = Blanciforti86[ set, ] )
print( estResultLAtrend )
summary( estResultLAtrend )

Blanciforti86$trend2 <- c( 0:( nrow( Blanciforti86 ) - 1 ) )^2
estResultLAtrend2 <- aidsEst( pNames, wNames, "xFood",
   shifterNames = c( "trend", "trend2" ), data = Blanciforti86[ set, ] )
print( estResultLAtrend2 )
summary( estResultLAtrend2 )


#####################################################
cat( paste( "\nRepeating the evaluation of different elasticity formulas",
   "of Green & Alston (1990): iterated AIDS\n" ) )
estResultAIDS <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86[ setWo1, ], method = "IL" )
print( estResultAIDS )
print( summary( estResultAIDS ) )
print( elas( estResultAIDS, method = "AIDS", quantNames = wNames ) )
# imposing restrictions via restrict.regMat
estResultAIDSTX <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86[ setWo1, ], method = "IL", restrict.regMat = TRUE )
print( estResultAIDSTX )
print( summary( estResultAIDSTX ) )
print( elas( estResultAIDSTX, method = "AIDS", quantNames = wNames ) )

## only homogeneity (no symmetry imposed)
estResultAIDShom <- aidsEst( pNames, wNames, "xFood", sym = FALSE,
   data = Blanciforti86[ setWo1, ], method = "IL" )
print( estResultAIDShom )
print( summary( estResultAIDShom ) )
print( elas( estResultAIDShom, method = "AIDS", quantNames = wNames ) )
# imposing restrictions via restrict.regMat
estResultAIDShomTX <- aidsEst( pNames, wNames, "xFood", sym = FALSE,
   data = Blanciforti86[ setWo1, ], method = "IL", restrict.regMat = TRUE )
print( estResultAIDShomTX )
print( summary( estResultAIDShomTX ) )
print( elas( estResultAIDShomTX, method = "AIDS", quantNames = wNames ) )

## unrestricted (no homogeneity and no symmetry imposed)
estResultAIDSunr <- aidsEst( pNames, wNames, "xFood", hom = FALSE, sym = FALSE,
   data = Blanciforti86[ setWo1, ], method = "IL" )
print( estResultAIDSunr )
print( summary( estResultAIDSunr ) )
print( elas( estResultAIDSunr, method = "AIDS", quantNames = wNames ) )
# imposing restrictions via restrict.regMat
estResultAIDSunrTX <- aidsEst( pNames, wNames, "xFood", hom = FALSE, sym = FALSE,
   data = Blanciforti86[ setWo1, ], method = "IL", restrict.regMat = TRUE )
print( estResultAIDSunrTX )
print( summary( estResultAIDSunrTX ) )
print( elas( estResultAIDSunrTX, method = "AIDS", quantNames = wNames ) )

## with NAs
estResultLaSNa <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86,
   priceIndex = "S" )
print( estResultLaSNa )
print( summary( estResultLaSNa ) )
print( elas( estResultLaSNa, method = "AIDS", quantNames = wNames ) )

estResultLaSlNa <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86,
   priceIndex = "SL" )
print( estResultLaSlNa )
print( summary( estResultLaSlNa ) )
print( elas( estResultLaSlNa, method = "AIDS", quantNames = wNames ) )

estResultLaLsNa <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86 )
print( estResultLaLsNa )
print( summary( estResultLaLsNa ) )
print( elas( estResultLaLsNa, method = "AIDS", quantNames = wNames ) )

estResultAIDSNa <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86, method = "IL" )
print( estResultAIDSNa )
print( summary( estResultAIDSNa ) )
print( elas( estResultAIDSNa, method = "AIDS", quantNames = wNames ) )


########## Elasticities ###############
cat( "\nAIDS: Elasticities\n" )
ela <- aidsElas( estResultAIDS$coef, wMeans, pMeans, method = "AIDS",
   coefVcov = estResultAIDS$coef$allcov, df = estResultAIDS$est$df )
print( ela )
print( summary( ela ) )
elaTX <- aidsElas( estResultAIDSTX$coef, wMeans, pMeans, method = "AIDS",
   coefVcov = estResultAIDSTX$coef$allcov, df = estResultAIDSTX$est$df )
print( elaTX )
print( summary( elaTX ) )
print( all.equal( ela, elaTX ) )

print( elas( estResultAIDS ) )
print( summary( elas( estResultAIDS ) ) )

print( elas( estResultAIDSTX ) )
print( summary( elas( estResultAIDSTX ) ) )


cat( "\nLA: Elasticity formula of non-linear AIDS\n" )
ela <- aidsElas( estResultLA$coef, wMeans, pMeans, method = "AIDS",
   coefVcov = estResultLA$coef$allcov, df = estResultLA$est$df )
print( ela )
print( summary( ela ) )
elaTX <- aidsElas( estResultLATX$coef, wMeans, pMeans, method = "AIDS",
   coefVcov = estResultLATX$coef$allcov, df = estResultLATX$est$df )
print( elaTX )
print( summary( elaTX ) )
print( all.equal( ela, elaTX ) )

print( elas( estResultLA, method = "AIDS" ) )
print( summary( elas( estResultLA, method = "AIDS" ) ) )

print( elas( estResultLATX, method = "AIDS" ) )
print( summary( elas( estResultLATX, method = "AIDS" ) ) )


cat( "\n********** Elasticities ***************" )
cat( "\nLA: Elasticity formula of Goddard or Chalfant\n" )
ela <- aidsElas( estResultLA$coef, wMeans, method = "Go",
   coefVcov = estResultLA$coef$allcov, df = estResultLA$est$df )
print( ela )
print( summary( ela ) )
ela <- aidsElas( estResultLA$coef, wMeans, method = "Ch",
   coefVcov = estResultLA$coef$allcov, df = estResultLA$est$df )
print( ela )
print( summary( ela ) )

print( elas( estResultLA, method = "Go" ) )
print( summary( elas( estResultLA ) ) )

print( elas( estResultLATX ) )
print( summary( elas( estResultLATX ) ) )


cat( "\nLA: Elasticity formula of Eales + Unnevehr\n" )
ela <- aidsElas( estResultLA$coef, wMeans, method = "EU" )
print( ela )

cat( "\nLA: Elasticity formula of Green + Alston\n" )
ela <- aidsElas( estResultLA$coef, wMeans, pMeans, method = "GA" )
print( ela )

cat( "\nLA: Elasticity formula of Buse\n" )
ela <- aidsElas( estResultLA$coef, wMeans, pMeans, method = "B1" )
print( ela )

cat( "\nLA: Elasticity formula of Buse (alternative formula)\n" )
ela <- aidsElas( estResultLA$coef, wMeans, pMeans, method = "B2" )
print( ela )


############# Price indices ##############
options( digits = 5 )
cat( "\n************** Price indices **************\n" )
cat( "\nStone index\n" )
pxS <- aidsPx( "S", pNames, shareNames = wNames, data = Blanciforti86 )
print( pxS )

cat( "\nStone index with lagged shares\n" )
pxSL <- aidsPx( "SL", pNames, shareNames = wNames, data = Blanciforti86 )
print( pxSL )

cat( "\nPaasche index\n" )
pxP <- aidsPx( "P", pNames, shareNames = wNames, data = Blanciforti86 )
print( pxP )

pxP2 <- aidsPx( "P", pNames, data = Blanciforti86, shareNames = wNames,
   base = row.names(Blanciforti86) == "1970" )
print( pxP2 )

cat( "\nLaspeyres index, simplified\n" )
pxLs <- aidsPx( "Ls", pNames, shareNames = wNames, data = Blanciforti86 )
print( pxLs )

pxLs2 <- aidsPx( "Ls", pNames, data = Blanciforti86, shareNames = wNames,
   base = c( 1:32 ) )
print( pxLs2 )

cat( "\nTornqvist index\n" )
pxT <- aidsPx( "T", pNames, shareNames = wNames, data = Blanciforti86 )
print( pxT )

pxT2 <- aidsPx( "T", pNames, data = Blanciforti86, shareNames = wNames,
   base = list( prices = rep( 100, 4 ), shares = rep( 0.25, 4 ) ) )
print( pxT2 )

cat( "\nTranslog index\n" )
pxTL <- aidsPx( "TL", pNames, shareNames = wNames, data = Blanciforti86,
   coef = c( list( alpha0 = 0 ), estResultLA$coef ) )
print( pxTL )

# Translog index with 1 demand shifter
pxTLtrend <- aidsPx( "TL", pNames, data = Blanciforti86,
   coef = c( list( alpha0 = 0 ), estResultLAtrend$coef ),
   shifterNames = c( "trend" ) )
print( pxTLtrend )

# Translog index with 2 demand shifters
pxTLtrend2 <- aidsPx( "TL", pNames, data = Blanciforti86,
   coef = c( list( alpha0 = 0 ), estResultLAtrend2$coef ),
   shifterNames = c( "trend", "trend2" ) )
print( pxTLtrend2 )


########### fitted values #################
options( digits = 3 )
fittedAIDS <- aidsCalc( pNames, "xFood", data = Blanciforti86[ -1, ],
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
fittedAIDSTX <- aidsCalc( pNames, "xFood", data = Blanciforti86[ -1, ],
   coef = estResultAIDSTX$coef )
print( fittedAIDSTX )
print( all.equal( fittedAIDS, fittedAIDSTX ) )

# LA-AIDS with Stone price price index with lagged shares
fittedLA <- aidsCalc( pNames, "xFood", data = Blanciforti86[ set, ],
   coef = estResultLA$coef, priceIndex = estResultLA$lnp )
print( fittedLA )
if( max( abs( fittedLA$shares[ -1, ] - estResultLA$wFitted[ setWo1, ] ),
      na.rm = TRUE ) > 1e-5 ) {
   stop( "fitted shares of LA-AIDS are wrong" )
}
if( max( abs( fittedLA$quant[ -1, ] - estResultLA$qFitted[ setWo1, ] ),
      na.rm = TRUE ) > 1e-5 ) {
   stop( "fitted quantities of LA-AIDS are wrong" )
}
fittedLATX <- aidsCalc( pNames, "xFood", data = Blanciforti86[ set, ],
   coef = estResultLATX$coef, priceIndex = estResultLATX$lnp )
print( fittedLATX )
print( all.equal( fittedLA, fittedLATX ) )

# LA-AIDS with Stone price index
# obsereved shares in the Stone price index
fittedLaSNa <- aidsCalc( pNames, "xFood", data = Blanciforti86[ set, ],
   coef = estResultLaSNa$coef, priceIndex = estResultLaSNa$lnp )
print( fittedLaSNa )
all.equal( fittedLaSNa$shares, estResultLaSNa$wFitted )
all.equal( fittedLaSNa$quant, estResultLaSNa$qFitted )
# fitted shares in the Stone price index
fittedLaSNa2 <- aidsCalc( pNames, "xFood", data = Blanciforti86[ set, ],
   coef = estResultLaSNa$coef, priceIndex = "S" )
print( fittedLaSNa2 )
B86LaSNa2 <- cbind( Blanciforti86[ set, c( pNames, "xFood" ) ],
   fittedLaSNa2$shares )
lnp <- aidsPx( "S", pNames, shareNames = c( "w1", "w2", "w3", "w4" ),
   data = B86LaSNa2 )
fittedLaSNa2b <- aidsCalc( pNames, "xFood", data = Blanciforti86[ set, ],
   coef = estResultLaSNa$coef, priceIndex = lnp )
all.equal( fittedLaSNa2, fittedLaSNa2b )

# LA-AIDS with simplified Laspeyres price index and NAs
fittedLaLsNa <- aidsCalc( pNames, "xFood", data = Blanciforti86[ set, ],
   coef = estResultLaLsNa$coef, priceIndex = estResultLaLsNa$lnp )
print( fittedLaLsNa )
all.equal( fittedLaLsNa$shares, estResultLaLsNa$wFitted[ set, ],
   check.attributes = FALSE )
all.equal( fittedLaLsNa$quant, estResultLaLsNa$qFitted[ set, ],
   check.attributes = FALSE )
fittedLaLsNa2 <- aidsCalc( pNames, "xFood", data = Blanciforti86[ set, ],
   coef = estResultLaLsNa$coef, priceIndex = "Ls",
   basePrices = as.numeric( Blanciforti86[ 1, pNames ] ),
   baseShares = as.numeric( Blanciforti86[ 1, wNames ] ) )
all.equal( estResultLaLsNa$wFitted, fittedLaLsNa2$shares,
   check.attributes = FALSE )
all.equal( estResultLaLsNa$qFitted, fittedLaLsNa2$quant,
   check.attributes = FALSE )


####### consistency ###################
# with observed expenditure shares
consist <- aidsTestConsist( pNames, "xFood", Blanciforti86[ set, ],
   coef = estResultAIDS$coef, shareNames = wNames )
print( consist )
class( consist ) <- NULL
print( consist )

# with fitted expenditure shares
consistFitted <- aidsTestConsist( pNames, totExpName = "xFood",
   data = Blanciforti86[ set, ], coef = estResultAIDS$coef )
print( consistFitted )
class( consistFitted ) <- NULL
print( consistFitted )


## replicating the LA-AIDS estimation of the SAS example
# loading data set
data( USMeatConsump )

# adding shifter variables for modeling seasonal effects
USMeatConsump$co1 <- cos( 1 / 2 * 3.14159 * USMeatConsump$t )
USMeatConsump$si1 <- sin( 1 / 2 * 3.14159 * USMeatConsump$t )

# Scaling prices by their means
USMeatConsump$beef_pm <- USMeatConsump$beef_p / mean( USMeatConsump$beef_p )
USMeatConsump$pork_pm <- USMeatConsump$pork_p / mean( USMeatConsump$pork_p )
USMeatConsump$chick_pm <- USMeatConsump$chick_p / mean( USMeatConsump$chick_p )
USMeatConsump$turkey_pm <- USMeatConsump$turkey_p / mean( USMeatConsump$turkey_p )

# Estimation of the model
meatModel <- aidsEst( c( "beef_pm", "pork_pm", "chick_pm", "turkey_pm" ),
   c( "beef_w", "pork_w", "chick_w", "turkey_w" ),
   "meat_exp", shifterNames = c( "co1", "si1", "t" ),
   method = "LA:S", data = USMeatConsump, maxiter=1000 )
meatModel
summary( meatModel )


## log likelihood values
logLik( estResultLA )
logLik( estResultLATX )
logLik( estResultLAhom )
logLik( estResultLAhomTX )
logLik( estResultLAunr )
logLik( estResultLAunrTX )
logLik( estResultLAtrend )
logLik( estResultLAtrend2 )
logLik( estResultAIDS )
logLik( estResultAIDSTX )
logLik( estResultAIDShom )
logLik( estResultAIDShomTX )
logLik( estResultAIDSunr )
logLik( estResultAIDSunrTX )
logLik( estResultLaSNa )
logLik( estResultLaSlNa )
logLik( estResultLaLsNa )
logLik( estResultAIDSNa )
logLik( meatModel )


## LR tests
lrtest( estResultLA, estResultLAhom, estResultLAunr, estResultLA )
lrtest( estResultLATX, estResultLAhomTX, estResultLAunrTX, estResultLATX )
lrtest( estResultLA, estResultLAtrend, estResultLAtrend2, estResultLA )
lrtest( estResultAIDSunr, estResultAIDShom, estResultAIDS, estResultAIDSunr )
lrtest( estResultAIDSunrTX, estResultAIDShomTX, estResultAIDSTX,
   estResultAIDSunrTX )


## comparing estimations results with different methods to impose restrictions
# estResultLA vs. estResultLATX
estResultLATX$call <- NULL
estResultLATX$est$call <- NULL
estResultLATX$est$restrict.regMat <- NULL
estResultLA$call <- NULL
estResultLA$est$call <- NULL
estResultLA$est$restrict.matrix <- NULL
estResultLA$est$restrict.rhs <- NULL
print( all.equal( estResultLA, estResultLATX ) )

# estResultLAhom vs. estResultLAhomTX
estResultLAhomTX$call <- NULL
estResultLAhomTX$est$call <- NULL
estResultLAhomTX$est$restrict.regMat <- NULL
estResultLAhom$call <- NULL
estResultLAhom$est$call <- NULL
estResultLAhom$est$restrict.matrix <- NULL
estResultLAhom$est$restrict.rhs <- NULL
print( all.equal( estResultLAhom, estResultLAhomTX ) )

# estResultLAunr vs. estResultLAunrTX
estResultLAunrTX$call <- NULL
estResultLAunrTX$est$call <- NULL
estResultLAunr$call <- NULL
estResultLAunr$est$call <- NULL
print( all.equal( estResultLAunr, estResultLAunrTX ) )

# estResultAIDS vs. estResultAIDSTX
estResultAIDSTX$call <- NULL
estResultAIDSTX$est$call <- NULL
estResultAIDSTX$est$restrict.regMat <- NULL
estResultAIDS$call <- NULL
estResultAIDS$est$call <- NULL
estResultAIDS$est$restrict.matrix <- NULL
estResultAIDS$est$restrict.rhs <- NULL
print( all.equal( estResultAIDS, estResultAIDSTX ) )

# estResultAIDShom vs. estResultAIDShomTX
estResultAIDShomTX$call <- NULL
estResultAIDShomTX$est$call <- NULL
estResultAIDShomTX$est$restrict.regMat <- NULL
estResultAIDShom$call <- NULL
estResultAIDShom$est$call <- NULL
estResultAIDShom$est$restrict.matrix <- NULL
estResultAIDShom$est$restrict.rhs <- NULL
print( all.equal( estResultAIDShom, estResultAIDShomTX ) )

# estResultAIDSunr vs. estResultAIDSunrTX
estResultAIDSunrTX$call <- NULL
estResultAIDSunrTX$est$call <- NULL
estResultAIDSunr$call <- NULL
estResultAIDSunr$est$call <- NULL
print( all.equal( estResultAIDSunr, estResultAIDSunrTX ) )
