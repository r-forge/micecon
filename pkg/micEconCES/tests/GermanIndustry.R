library( "micEconCES" )
data( "GermanIndustry" )

print( GermanIndustry )

GermanIndustry$time <- GermanIndustry$year - 1963

xNames <-  c( "K", "E", "A" )

b <- c( "gamma" = 38, "lambda" = 0.0222, "delta_1" = 0.5, "delta_2" = 0.1, 
   "rho_1" = 0.5300, "rho" = 0.1813 )

GermanIndustry$YCalc <- cesCalc( xNames = xNames, tName = "time",
   data = GermanIndustry, coef = b, nested = TRUE )

GermanIndustry$YCalc 


################# econometric estimation with cesEst ##########################

## Nelder-Mead
cesNm <- cesEst( "Y", xNames, tName = "time", data = GermanIndustry,
   method = "NM", control = list( maxit = 5000 ) )
print.default( cesNm )
print( cesNm )
summary( cesNm )

## Conjugate Gradients
cesCg <- cesEst( "Y", xNames, tName = "time", data = GermanIndustry,
   method = "CG", control = list( maxit = 1000 ) )
print.default( cesCg )
print( cesCg )
summary( cesCg )

## Simulated Annealing
cesSann <- cesEst( "Y", xNames, tName = "time", data = GermanIndustry,
   method = "SANN", control = list( maxit = 20000 ) )
print.default( cesSann )
print( cesSann )
summary( cesSann )

## BFGS
cesBfgs <- cesEst( "Y", xNames, tName = "time", data = GermanIndustry,
   method = "BFGS", control = list( maxit = 5000 ) )
print.default( cesBfgs )
print( cesBfgs )
summary( cesBfgs )

## L-BFGS-B
cesBfgsCon <- cesEst( "Y", xNames, tName = "time", data = GermanIndustry,
   method = "L-BFGS-B", control = list( maxit = 5000 )  )
print.default( cesBfgsCon )
print( cesBfgsCon )
summary( cesBfgsCon )

## Levenberg-Marquardt
cesLm <- cesEst( "Y", xNames, tName = "time", data = GermanIndustry,
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
print.default( cesLm )
print( cesLm )
summary( cesLm )

## Newton-type
cesNewton <- cesEst( "Y", xNames, tName = "time", data = GermanIndustry,
   method = "Newton", iterlim = 500 )
print.default( cesNewton )
print( cesNewton )
summary( cesNewton )

## PORT
cesPort <- cesEst( "Y", xNames, tName = "time", data = GermanIndustry,
   method = "PORT", control = list( eval.max = 1000, iter.max = 1000 ) )
print.default( cesPort )
print( cesPort )
summary( cesPort )

## DE
cesDe <- cesEst( "Y", xNames, tName = "time", data = GermanIndustry,
   method = "DE", control = DEoptim.control( trace = FALSE, NP = 90 ) )
print.default( cesDe )
print( cesDe )
summary( cesDe )

## nls
try( cesNls <- cesEst( "Y", xNames, tName = "time", data = GermanIndustry,
   vrs = TRUE, method = "nls" ) )


################ estimation with rho_1 and/or rho fixed #######################

## Nelder-Mead, rho_1 and rho fixed
cesNmRho2 <- cesEst( "Y", xNames, tName = "time", data = GermanIndustry,
   method = "NM", rho1 = 0.1813, rho = 0.5300, 
   control = list( maxit = 5000 ) )
print.default( cesNmRho2 )
print( cesNmRho2 )
summary( cesNmRho2 )

## BFGS, rho_1 and rho fixed
cesBfgsRho2 <- cesEst( "Y", xNames, tName = "time", data = GermanIndustry,
   method = "BFGS", rho1 = 0.1813, rho = 0.5300, 
   control = list( maxit = 5000 ) )
print.default( cesBfgsRho2 )
print( cesBfgsRho2 )
summary( cesBfgsRho2 )

## L-BFGS-B, rho_1 and rho fixed
cesBfgsConRho2 <- cesEst( "Y", xNames, tName = "time", data = GermanIndustry,
   method = "L-BFGS-B", rho1 = 0.1813, rho = 0.5300, 
   control = list( maxit = 5000 )  )
print.default( cesBfgsConRho2 )
print( cesBfgsConRho2 )
summary( cesBfgsConRho2 )

## Levenberg-Marquardt, rho fixed
cesLmRho <- cesEst( "Y", xNames, tName = "time", data = GermanIndustry, 
   rho = 0.5300, control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
print.default( cesLmRho )
print( cesLmRho )
summary( cesLmRho )

## Levenberg-Marquardt, rho_1 fixed
cesLmRho1 <- cesEst( "Y", xNames, tName = "time", data = GermanIndustry, 
   rho1 = 0.1813, control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
print.default( cesLmRho1 )
print( cesLmRho1 )
summary( cesLmRho1 )

## Levenberg-Marquardt, rho_1 and rho fixed
cesLmRho2 <- cesEst( "Y", xNames, tName = "time", data = GermanIndustry, 
   rho1 = 0.1813, rho = 0.5300,
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
print.default( cesLmRho2 )
print( cesLmRho2 )
summary( cesLmRho2 )

## Newton-type, rho_1 and rho fixed
cesNewtonRho2 <- cesEst( "Y", xNames, tName = "time", data = GermanIndustry,
   method = "Newton", rho1 = 0.1813, rho = 0.5300, iterlim = 500 )
print.default( cesNewtonRho2 )
print( cesNewtonRho2 )
summary( cesNewtonRho2 )

## PORT, rho_1 and rho fixed
cesPortRho2 <- cesEst( "Y", xNames, tName = "time", data = GermanIndustry,
   method = "PORT", rho1 = 0.1813, rho = 0.5300, 
   control = list( eval.max = 1000, iter.max = 1000 ) )
print.default( cesPortRho2 )
print( cesPortRho2 )
summary( cesPortRho2 )

## DE, rho_1 and rho fixed
cesDeRho2 <- cesEst( "Y", xNames, tName = "time", data = GermanIndustry,
   method = "DE", rho1 = 0.1813, rho = 0.5300, 
   control = DEoptim.control( trace = FALSE, NP = 90 ) )
print.default( cesDeRho2 )
print( cesDeRho2 )
summary( cesDeRho2 )


