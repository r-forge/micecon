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

# estimation
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

# estimation with Epanechnikov kernel
npseed( 123 )
estResultEpa <- npregHom( "qVarInput",
   xNames = c( "pOutput", "pVarInput", "pLabor", "land" ),
   data = germanFarms, homWeights = weights, ckertype="epanechnikov" )
print( estResultEpa )
all.equal( estResultEpa$grad[ , "pOutput" ] * germanFarms$pOutput +
   estResultEpa$grad[ , "pVarInput" ] * germanFarms$pVarInput,
   - estResultEpa$grad[ , "pLabor" ] * germanFarms$pLabor )

# estimation with manual bandwidth selection
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
