library( "miscTools" )
library( "graphicsQC" )

set.seed( 123 )
x <- 1000 * runif( 25 )
y <- 10 + 3 * x + 100 * rnorm( 25 )
ols <- lm( y ~ x )

obsLab <- "observed"
fitLab <- "fitted"
xyString <- "xy"

expressions <- c( "compPlot( y, fitted( ols ) )",
   "compPlot( y, fitted( ols ), lim = c( 0, 3500 ) )",
   "compPlot( y, fitted( ols ), pch = 20 )",
   "compPlot( y, fitted( ols ), xlab = obsLab, ylab = fitLab )",
   "compPlot( y, fitted( ols ), log = xyString )" )

# generate control plots:                  
# plotExpr( expressions, path = "tests/controlPlots", prefix = "compPlot", clear = TRUE )

plotExpr( expressions, path = "testPlots", prefix = "compPlot" )

compare( "testPlots", "controlPlots" )

