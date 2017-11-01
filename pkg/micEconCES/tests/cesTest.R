library( "micEconCES" )
options( digits = 3 )

# load data
data( germanFarms, package = "micEcon" )
# output quantity:
germanFarms$qOutput <- germanFarms$vOutput / germanFarms$pOutput
# quantity of intermediate inputs
germanFarms$qVarInput <- germanFarms$vVarInput / germanFarms$pVarInput


## CES: Land & Labor (Nelder-Mead)
cesLandLabor <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms,
   method = "NM", returnGrad = TRUE )
print.default( cesLandLabor ) 
print( cesLandLabor )
summary( cesLandLabor )
summary( cesLandLabor, ela = FALSE )
print( summary( cesLandLabor ), ela = FALSE )
summary( cesLandLabor )$elaCov
coef( cesLandLabor ) 
vcov( cesLandLabor ) 
coef( summary( cesLandLabor ) )
fitted( cesLandLabor )
residuals( cesLandLabor )

# variable returns to scale (Nelder-Mead)
cesLandLaborVrs <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms,
   vrs = TRUE, method = "Nelder-Mead", returnGrad = TRUE )
print.default( cesLandLaborVrs )
print( cesLandLaborVrs )
summary( cesLandLaborVrs )
coef( cesLandLaborVrs )
vcov( cesLandLaborVrs )
coef( summary( cesLandLaborVrs ) )
fitted( cesLandLaborVrs )
residuals( cesLandLaborVrs )

# using the CG optimization method
cesLandLaborCg <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms,
   method = "CG" )
print.default( cesLandLaborCg )
print( cesLandLaborCg )
summary( cesLandLaborCg )
coef( cesLandLaborCg )
vcov( cesLandLaborCg )
coef( summary( cesLandLaborCg ) )
fitted( cesLandLaborCg )
residuals( cesLandLaborCg )

# using the SANN optimization method
cesLandLaborSann <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms,
   method = "SANN", returnGrad = TRUE )
print.default( cesLandLaborSann )
print( cesLandLaborSann )
summary( cesLandLaborSann )
coef( cesLandLaborSann )
vcov( cesLandLaborSann )
coef( summary( cesLandLaborSann ) )
fitted( cesLandLaborSann )
residuals( cesLandLaborSann )

# using the BFGS optimization method
cesLandLaborBfgs <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms,
   method = "BFGS" )
print.default( cesLandLaborBfgs )
print( cesLandLaborBfgs )
summary( cesLandLaborBfgs )
coef( cesLandLaborBfgs )
vcov( cesLandLaborBfgs )
coef( summary( cesLandLaborBfgs ) )
fitted( cesLandLaborBfgs )
residuals( cesLandLaborBfgs )

# using the L-BFGS-B optimization method with constrained parameters
cesLandLaborBfgsCon <- cesEst( "qOutput", c( "land", "qLabor" ),
   germanFarms, method = "L-BFGS-B", returnGrad = TRUE )
print.default( cesLandLaborBfgsCon )
print( cesLandLaborBfgsCon )
summary( cesLandLaborBfgsCon )
coef( cesLandLaborBfgsCon )
vcov( cesLandLaborBfgsCon )
coef( summary( cesLandLaborBfgsCon ) )
fitted( cesLandLaborBfgsCon )
residuals( cesLandLaborBfgsCon )

# Kmenta approximation with CRS
cesLandLaborKmentaCrs <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms,
   method = "Kmenta" )
print.default( cesLandLaborKmentaCrs )
print( cesLandLaborKmentaCrs )
summary( cesLandLaborKmentaCrs )
coef( cesLandLaborKmentaCrs )
vcov( cesLandLaborKmentaCrs )
coef( summary( cesLandLaborKmentaCrs ) )
fitted( cesLandLaborKmentaCrs )
residuals( cesLandLaborKmentaCrs )

# Kmenta approximation with VRS
cesLandLaborKmenta <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms,
   vrs = TRUE, method = "Kmenta", returnGrad = TRUE )
print.default( cesLandLaborKmenta )
print( cesLandLaborKmenta )
summary( cesLandLaborKmenta )
coef( cesLandLaborKmenta )
vcov( cesLandLaborKmenta )
coef( summary( cesLandLaborKmenta ) )
fitted( cesLandLaborKmenta )
residuals( cesLandLaborKmenta )

# using the Levenberg-Marquardt optimization method
cesLandLaborLm <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms,
   method = "LM", control = nls.lm.control( maxiter = 200 ) )
print.default( cesLandLaborLm )
print( cesLandLaborLm )
summary( cesLandLaborLm )
coef( cesLandLaborLm )
vcov( cesLandLaborLm )
coef( summary( cesLandLaborLm ) )
fitted( cesLandLaborLm )
residuals( cesLandLaborLm )

# using the Newton-type optimization method implemented in nlm()
cesLandLaborNewton <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms,
   method = "Newton", returnGrad = TRUE )
print.default( cesLandLaborNewton )
print( cesLandLaborNewton )
summary( cesLandLaborNewton )
coef( cesLandLaborNewton )
vcov( cesLandLaborNewton )
coef( summary( cesLandLaborNewton ) )
fitted( cesLandLaborNewton )
residuals( cesLandLaborNewton )

# using the PORT optimization rountine implemented in nlminb(), UNconstrained
cesLandLaborPort <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms,
   method = "PORT", lower = -Inf, upper = Inf )
print.default( cesLandLaborPort )
print( cesLandLaborPort )
summary( cesLandLaborPort )
summary( cesLandLaborPort, ela = FALSE )
print( summary( cesLandLaborPort ), ela = FALSE )
summary( cesLandLaborPort )$elaCov
coef( cesLandLaborPort )
vcov( cesLandLaborPort )
coef( summary( cesLandLaborPort ) )
fitted( cesLandLaborPort )
residuals( cesLandLaborPort )

# using the PORT optimization rountine implemented in nlminb(), constrained
cesLandLaborPortCon <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms,
   method = "PORT", returnGrad = TRUE )
print.default( cesLandLaborPortCon )
print( cesLandLaborPortCon )
summary( cesLandLaborPortCon )
coef( cesLandLaborPortCon )
vcov( cesLandLaborPortCon )
coef( summary( cesLandLaborPortCon ) )
fitted( cesLandLaborPortCon )
residuals( cesLandLaborPortCon )

# using the PORT optimization rountine implemented in nlminb(), constrained by hand
cesLandLaborPortCon2 <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms,
   method = "PORT", lower = c( 0.1, 0.1, -0.5 ), upper = c( 300, 0.9, 4 ) )
print.default( cesLandLaborPortCon2 )
print( cesLandLaborPortCon2 )
summary( cesLandLaborPortCon2 )
coef( cesLandLaborPortCon2 )
vcov( cesLandLaborPortCon2 )
coef( summary( cesLandLaborPortCon2 ) )
fitted( cesLandLaborPortCon2 )
residuals( cesLandLaborPortCon2 )

# using the PORT optimization rountine implemented in nlminb(), VRS, UNconstrained
cesLandLaborPortVrs <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms,
   method = "PORT", vrs = TRUE, lower = -Inf, upper = Inf )
print.default( cesLandLaborPortVrs )
print( cesLandLaborPortVrs )
summary( cesLandLaborPortVrs )
coef( cesLandLaborPortVrs )
vcov( cesLandLaborPortVrs )
coef( summary( cesLandLaborPortVrs ) )
fitted( cesLandLaborPortVrs )
residuals( cesLandLaborPortVrs )

# using the PORT optimization rountine implemented in nlminb(), VRS, constrained
cesLandLaborPortVrsCon <- cesEst( "qOutput", c( "land", "qLabor" ),
   germanFarms, method = "PORT", vrs = TRUE )
print.default( cesLandLaborPortVrsCon )
print( cesLandLaborPortVrsCon )
summary( cesLandLaborPortVrsCon )
coef( cesLandLaborPortVrsCon )
vcov( cesLandLaborPortVrsCon )
coef( summary( cesLandLaborPortVrsCon ) )
fitted( cesLandLaborPortVrsCon )
residuals( cesLandLaborPortVrsCon )

# using the DE optimization method implemented in DEoptim(), CRS
cesLandLaborDe <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms,
   method = "DE", control = DEoptim.control( trace = FALSE ), 
   returnGrad = TRUE )
print.default( cesLandLaborDe )
print( cesLandLaborDe )
summary( cesLandLaborDe )
coef( cesLandLaborDe )
vcov( cesLandLaborDe )
coef( summary( cesLandLaborDe ) )
fitted( cesLandLaborDe )
residuals( cesLandLaborDe )

# using the DE optimization method implemented in DEoptim(), VRS
cesLandLaborDeVrs <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms,
   method = "DE", vrs = TRUE, control = DEoptim.control( trace = FALSE ) )
print.default( cesLandLaborDeVrs )
print( cesLandLaborDeVrs )
summary( cesLandLaborDeVrs )
coef( cesLandLaborDeVrs )
vcov( cesLandLaborDeVrs )
coef( summary( cesLandLaborDeVrs ) )
fitted( cesLandLaborDeVrs )
residuals( cesLandLaborDeVrs )

# using the DE optimization method implemented in DEoptim(), CRS, user-specified bounds
cesLandLaborDe2 <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms,
   method = "DE", control = DEoptim.control( trace = FALSE ),
   lower = c( 1, 0.01, -0.99 ), upper = c( 50, 0.99, 3 ) )
print.default( cesLandLaborDe2 )
print( cesLandLaborDe2 )
summary( cesLandLaborDe2 )
coef( cesLandLaborDe2 )
vcov( cesLandLaborDe2 )
coef( summary( cesLandLaborDe2 ) )
fitted( cesLandLaborDe2 )
residuals( cesLandLaborDe2 )

# Land & Labor with constant returns to scale (nls)
try( cesLandLaborNls <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms,
   method = "nls" ) )

# Land & Labor with variable returns to scale (nls)
try( cesLandLaborNls <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms,
   method = "nls", vrs = TRUE ) )


## CES: Land & Intermediate Inputs (Nelder-Mead)
cesLandInt <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms,
   method = "NM" )
print.default( cesLandInt )
print( cesLandInt )
summary( cesLandInt )
coef( summary( cesLandInt ) )

# variable returns to scale (Nelder-Mead)
cesLandIntVrs <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms,
   vrs = TRUE, method = "Nelder-Mead" )
print.default( cesLandIntVrs )
print( cesLandIntVrs )
summary( cesLandIntVrs )
coef( summary( cesLandIntVrs ) )

# using the CG optimization method
cesLandIntCg <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms,
   method = "CG" )
print.default( cesLandIntCg )
print( cesLandIntCg )
summary( cesLandIntCg )
coef( summary( cesLandIntCg ) )

# using the SANN optimization method
cesLandIntSann <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms,
   method = "SANN", random.seed = 234 )
print.default( cesLandIntSann )
print( cesLandIntSann )
summary( cesLandIntSann )
coef( summary( cesLandIntSann ) )

# using the BFGS optimization method
cesLandIntBfgs <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms,
   method = "BFGS" )
print.default( cesLandIntBfgs )
print( cesLandIntBfgs )
summary( cesLandIntBfgs )
coef( summary( cesLandIntBfgs ) )

# using the L-BFGS-B optimization method with constrained parameters
cesLandIntBfgsCon <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms,
   method = "L-BFGS-B" )
print.default( cesLandIntBfgsCon )
print( cesLandIntBfgsCon )
summary( cesLandIntBfgsCon )
coef( summary( cesLandIntBfgsCon ) )

# Kmenta approximation with CRS
cesLandIntKmentaCrs <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms,
   method = "Kmenta" )
print.default( cesLandIntKmentaCrs )
print( cesLandIntKmentaCrs )
summary( cesLandIntKmentaCrs )
coef( summary( cesLandIntKmentaCrs ) )

# Kmenta approximation with VRS
cesLandIntKmenta <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms,
   method = "Kmenta", vrs = TRUE )
print.default( cesLandIntKmenta )
print( cesLandIntKmenta )
summary( cesLandIntKmenta )
coef( summary( cesLandIntKmenta ) )

# using the Levenberg-Marquardt optimization method
cesLandIntLm <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms,
   method = "LM", control = nls.lm.control( maxiter = 200 ) )
print.default( cesLandIntLm )
print( cesLandIntLm )
summary( cesLandIntLm )
coef( summary( cesLandIntLm ) )

# using the Newton-type optimization method implemented in nlm()
cesLandIntNewton <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms,
   method = "Newton" )
print.default( cesLandIntNewton )
print( cesLandIntNewton )
summary( cesLandIntNewton )
coef( summary( cesLandIntNewton ) )

# using the PORT optimization rountine implemented in nlminb(), UNconstrained
cesLandIntPort <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms,
   method = "PORT", lower = -Inf, upper = Inf )
print.default( cesLandIntPort )
print( cesLandIntPort )
summary( cesLandIntPort )
coef( summary( cesLandIntPort ) )

# using the PORT optimization rountine implemented in nlminb(), constrained
cesLandIntPortCon <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms,
   method = "PORT" )
print.default( cesLandIntPortCon )
print( cesLandIntPortCon )
summary( cesLandIntPortCon )
coef( summary( cesLandIntPortCon ) )

# using the PORT optimization rountine implemented in nlminb(), VRS, UNconstrained
cesLandIntPortVrs <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms,
   method = "PORT", vrs = TRUE, lower = -Inf, upper = Inf )
print.default( cesLandIntPortVrs )
print( cesLandIntPortVrs )
summary( cesLandIntPortVrs )
coef( summary( cesLandIntPortVrs ) )

# using the PORT optimization rountine implemented in nlminb(), VRS, constrained
cesLandIntPortVrsCon <- cesEst( "qOutput", c( "land", "qVarInput" ),
   germanFarms, method = "PORT", vrs = TRUE )
print.default( cesLandIntPortVrsCon )
print( cesLandIntPortVrsCon )
summary( cesLandIntPortVrsCon )
coef( summary( cesLandIntPortVrsCon ) )

# using the DE optimization method implemented in DEoptim()
cesLandIntDe <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms,
   method = "DE", control = DEoptim.control( trace = FALSE ) )
print.default( cesLandIntDe )
print( cesLandIntDe )
summary( cesLandIntDe )
coef( summary( cesLandIntDe ) )

# using the DE optimization method implemented in DEoptim(), VRS
cesLandIntDeVrs <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms,
   method = "DE", vrs = TRUE, control = DEoptim.control( trace = FALSE ) )
print.default( cesLandIntDeVrs )
print( cesLandIntDeVrs )
summary( cesLandIntDeVrs )
coef( summary( cesLandIntDeVrs ) )

# using the DE optimization method implemented in DEoptim(), CRS, user-specified bounds
cesLandLaborDe2 <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms,
   method = "DE", control = DEoptim.control( trace = FALSE ),
   lower = c( 1, 0.01, -0.99 ), upper = c( 50, 0.99, 3 ) )
print.default( cesLandLaborDe2 )
print( cesLandLaborDe2 )
summary( cesLandLaborDe2 )
coef( cesLandLaborDe2 )
vcov( cesLandLaborDe2 )
coef( summary( cesLandLaborDe2 ) )
fitted( cesLandLaborDe2 )
residuals( cesLandLaborDe2 )

# constant returns to scale (nls)
try( cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms,
   method = "nls" ) )

# variable returns to scale (nls)
try( cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms,
   vrs = TRUE, method = "nls" ) )


############  cesCalc  ################
outLandLabor <- cesCalc( c( "land", "qLabor" ), germanFarms,
   coef( cesLandLabor ) )
print( outLandLabor )
all.equal( outLandLabor, cesCalc( c( "land", "qLabor" ), germanFarms,
   coef( cesLandLabor )[ c( 2, 3, 1 ) ] ) )
all.equal( outLandLabor, cesCalc( c( "land", "qLabor" ), germanFarms,
   unname( coef( cesLandLabor ) ) ) )

outLandLaborVrs <- cesCalc( c( "land", "qLabor" ), germanFarms,
   coef( cesLandLaborVrs ) )
print( outLandLaborVrs )
all.equal( outLandLaborVrs, cesCalc( c( "land", "qLabor" ), germanFarms,
   coef( cesLandLaborVrs )[ c( 3, 1, 4, 2 ) ] ) )
all.equal( outLandLaborVrs, cesCalc( c( "land", "qLabor" ), germanFarms,
   unname( coef( cesLandLaborVrs ) ) ) )
