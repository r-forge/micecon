
data( germanFarms )

germanFarms$qOutput   <- germanFarms$vOutput / germanFarms$pOutput
germanFarms$qVarInput <- -germanFarms$vVarInput / germanFarms$pVarInput
germanFarms$qLabor    <- -germanFarms$qLabor
germanFarms$time      <- c( 0:19 )

pNamesT <- c( "pOutput", "pVarInput", "pLabor" )
qNamesT <- c( "qOutput", "qVarInput", "qLabor" )
fNamesT <- c( "land", "time" )

estResult <- snqProfitEst( pNamesT, qNamesT, "land", data = germanFarms )
print( estResult )

##############################################
estResult <- snqProfitEst( pNamesT, qNamesT, fNamesT, data=germanFarms )
print( estResult )

estResultCalc <- snqProfitCalc( pNamesT, fNamesT, estResult$estData,
   estResult$weights, estResult$coef )
print( estResultCalc )
if( max( abs( estResultCalc - estResult$fitted ) ) > 1e-5 ) {
   stop( "values from snqProfitCalc are not equal to fitted values." )
}

estResultEla <- snqProfitEla( estResult$coef$beta,
   estResult$estData[ 20, pNames ], estResult$estData[ 20, qNames ],
   estResult$weights )
print( estResultEla )

estResultHderiv <- snqProfitHderiv( estResult$pMean, estResult$weights, 2 )
print( estResultHderiv )

estResultHessian <- snqProfitHessian( estResult$coef$beta,
   estResult$estData[ 20, pNames ], estResult$weights )
print( estResultHessian )

estResultShadowprices <- snqProfitShadowPrices( pNamesT, fNamesT, estResult$estData,
   estResult$weights, estResult$coef )
print( estResultShadowprices )

####################################################
estResult <- snqProfitEst( pNamesT, qNamesT, fNamesT, data=germanFarms, form = 1 )
print( estResult )

estResultCalc <- snqProfitCalc( pNamesT, fNamesT, estResult$estData,
   estResult$weights, estResult$coef, form = 1 )
print( estResultCalc )
if( max( abs( estResultCalc - estResult$fitted ) ) > 1e-5 ) {
   #stop( "values from snqProfitCalc are not equal to fitted values." )
}

estResultEla <- snqProfitEla( estResult$coef$beta,
   estResult$estData[ 20, pNames ], estResult$estData[ 20, qNames ],
   estResult$weights )
print( estResultEla )

estResultHderiv <- snqProfitHderiv( estResult$pMean, estResult$weights, 2 )
print( estResultHderiv )

estResultHessian <- snqProfitHessian( estResult$coef$beta,
   estResult$estData[ 20, pNames ], estResult$weights )
print( estResultHessian )

estResultShadowprices <- snqProfitShadowPrices( pNamesT, fNamesT, estResult$estData,
   estResult$weights, estResult$coef, form = 1 )
print( estResultShadowprices )
