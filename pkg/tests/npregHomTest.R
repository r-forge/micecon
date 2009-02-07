# load micEcon package
library( micEcon )

# load data
data( germanFarms )
# output quantity
germanFarms$qOutput <- germanFarms$vOutput / germanFarms$pOutput
# quantity of variable inputs
germanFarms$qVarInput <- germanFarms$vVarInput / germanFarms$pVarInput
# a time trend to account for technical progress:
germanFarms$time <- c(1:20)

# weights to impose normalize prices
weights <- c(
   pOutput = mean( germanFarms$qOutput ),
   pVarInput = mean( germanFarms$qVarInput ),
   pLabor = mean( germanFarms$qLabor ) )
weights <- weights / sum( weights )

# estimation (restricted gradients)
npseed( 123 )
estResult <- npregHom( "qVarInput",
   xNames = c( "pOutput", "pVarInput", "pLabor", "land" ),
   data = germanFarms, homWeights = weights )
print( estResult )
all.equal( estResult$grad[ , "pOutput" ] * germanFarms$pOutput +
   estResult$grad[ , "pVarInput" ] * germanFarms$pVarInput,
   - estResult$grad[ , "pLabor" ] * germanFarms$pLabor )
# different normalized variable omitted
npseed( 123 )
estResult2 <- npregHom( "qVarInput",
   xNames = c( "pOutput", "pVarInput", "pLabor", "land" ),
   data = germanFarms, homWeights = weights[ c( 3, 2, 1 ) ] )
print( estResult2 )
all.equal( estResult2$grad[ , "pOutput" ] * germanFarms$pOutput +
   estResult2$grad[ , "pVarInput" ] * germanFarms$pVarInput,
   - estResult2$grad[ , "pLabor" ] * germanFarms$pLabor )

# estimation (gradients not restricted)
npseed( 123 )
estResultAll <- npregHom( "qVarInput",
   xNames = c( "pOutput", "pVarInput", "pLabor", "land" ),
   data = germanFarms, homWeights = weights, restrictGrad = FALSE )
print( estResultAll )
all.equal( estResultAll$grad[ , "pOutput" ] * germanFarms$pOutput +
   estResultAll$grad[ , "pVarInput" ] * germanFarms$pVarInput,
   - estResultAll$grad[ , "pLabor" ] * germanFarms$pLabor )
# different order of weights
npseed( 123 )
estResultAll2 <- npregHom( "qVarInput",
   xNames = c( "pOutput", "pVarInput", "pLabor", "land" ),
   data = germanFarms, homWeights = weights[ c( 3, 2, 1 ) ],
   restrictGrad = FALSE )
all.equal( estResultAll$grad, estResultAll2$grad, tolerance = 1e-6 )


# estimation with Epanechnikov kernel (restricted gradients)
npseed( 123 )
estResultEpa <- npregHom( "qVarInput",
   xNames = c( "pOutput", "pVarInput", "pLabor", "land" ),
   data = germanFarms, homWeights = weights, ckertype="epanechnikov" )
print( estResultEpa )
all.equal( estResultEpa$grad[ , "pOutput" ] * germanFarms$pOutput +
   estResultEpa$grad[ , "pVarInput" ] * germanFarms$pVarInput,
   - estResultEpa$grad[ , "pLabor" ] * germanFarms$pLabor )

# estimation with Epanechnikov kernel (gradients not restricted)
npseed( 123 )
estResultEpaAll <- npregHom( "qVarInput",
   xNames = c( "pOutput", "pVarInput", "pLabor", "land" ),
   data = germanFarms, homWeights = weights, restrictGrad = FALSE,
   ckertype="epanechnikov" )
print( estResultEpaAll )
all.equal( estResultEpaAll$grad[ , "pOutput" ] * germanFarms$pOutput +
   estResultEpaAll$grad[ , "pVarInput" ] * germanFarms$pVarInput,
   - estResultEpaAll$grad[ , "pLabor" ] * germanFarms$pLabor )
# different order of weights
npseed( 123 )
estResultEpaAll2 <- npregHom( "qVarInput",
   xNames = c( "pOutput", "pVarInput", "pLabor", "land" ),
   data = germanFarms, homWeights = weights[ c( 3, 2, 1 ) ],
   restrictGrad = FALSE, ckertype="epanechnikov" )
all.equal( estResultEpaAll$grad, estResultEpaAll2$grad, tolerance = 1e-6 )


# estimation with manual bandwidth selection (restricted gradients)
estResultMan <- npregHom( "qVarInput",
   xNames = c( "pOutput", "pVarInput", "pLabor", "land" ),
   data = germanFarms, homWeights = weights, bws = rep( 1, 3 ),
   bwscaling = TRUE )
print( estResultMan )
all.equal( estResultMan$grad[ , "pOutput" ] * germanFarms$pOutput +
   estResultMan$grad[ , "pVarInput" ] * germanFarms$pVarInput,
   - estResultMan$grad[ , "pLabor" ] * germanFarms$pLabor )
# different normalized variable omitted
estResultMan2 <- npregHom( "qVarInput",
   xNames = c( "pOutput", "pVarInput", "pLabor", "land" ),
   data = germanFarms, homWeights = weights[ c( 3, 2, 1 ) ],
   bws = rep( 1, 3 ), bwscaling = TRUE )
print( estResultMan2 )
all.equal( estResultMan2$grad[ , "pOutput" ] * germanFarms$pOutput +
   estResultMan2$grad[ , "pVarInput" ] * germanFarms$pVarInput,
   - estResultMan2$grad[ , "pLabor" ] * germanFarms$pLabor )

# estimation with manual bandwidth selection (gradients not restricted)
estResultManAll <- npregHom( "qVarInput",
   xNames = c( "pOutput", "pVarInput", "pLabor", "land" ),
   data = germanFarms, homWeights = weights, restrictGrad = FALSE,
   bws = rep( 1, 4 ), bwscaling = TRUE )
print( estResultManAll )
all.equal( estResultManAll$grad[ , "pOutput" ] * germanFarms$pOutput +
   estResultManAll$grad[ , "pVarInput" ] * germanFarms$pVarInput,
   - estResultManAll$grad[ , "pLabor" ] * germanFarms$pLabor )
# different order of weights
estResultManAll2 <- npregHom( "qVarInput",
   xNames = c( "pOutput", "pVarInput", "pLabor", "land" ),
   data = germanFarms, homWeights = weights[ c( 3, 2, 1 ) ],
   restrictGrad = FALSE, bws = rep( 1, 4 ), bwscaling = TRUE )
all.equal( estResultManAll$grad, estResultManAll2$grad )
