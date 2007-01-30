heckit <- function( selection, formula,
                   data=sys.frame(sys.parent()),
                   inst = NULL,
                   print.level = 0 ) {
   # What is the role of na.action here?  We cannot use na.omit -- we must not omit the observation
   # where outcome is not observed.  na-s cannot be passed either.
   # However, we can (and should?) omit the na-s in explanatory and probit outcomes.  This needs
   # a bit of refinement.
   if( class( formula ) != "formula" ) {
      stop( "argument 'formula' must be a formula" )
   } else if( length( formula ) != 3 ) {
      stop( "argument 'formula' must be a 2-sided formula" )
   } else if( "invMillsRatio" %in% all.vars( formula ) ) {
      stop( "argument 'formula' may not include a variable name",
         " 'invMillsRatio'" )
   } else if( class( selection ) != "formula" ) {
      stop( "argument 'selection' must be a formula" )
   } else if( length( selection ) != 3 ) {
      stop( "argument 'selection' must be a 2-sided formula" )
   } else if( "probit" %in% substr( all.vars( selection ), 1, 6 ) ) {
      stop( "argument 'selection' may not include a variable",
         " names starting with 'probit'" )
   } else if( !is.null( inst ) ) {
      if( class ( inst ) != "formula" || length( inst ) != 2 ) {
         stop( "argument 'inst' must be a 1-sided formula" )
      }
   }

   result <- list()
   ## Now extract model frames etc.
   ## Selection equation
   mf <- match.call(expand.dots = FALSE)
   m <- match(c("selection", "data", "subset", "weights", "na.action",
                "offset"), names(mf), 0)
   mfS <- mf[c(1, m)]
   mfS$drop.unused.levels <- TRUE
   mfS[[1]] <- as.name("model.frame")
   names(mfS)[2] <- "formula"
                                        # model.frame requires the parameter to
                                        # be 'formula'
   mfS <- eval(mfS, parent.frame())
   mtS <- attr(mfS, "terms")
   firstStepData <- model.matrix(mtS, mfS)
   NXS <- ncol(firstStepData)
   probitEndogenous <- factor(model.response(mfS, "numeric"))
   probitLevels <- levels( as.factor( probitEndogenous ) )
   if( length( probitLevels ) != 2 ) {
      stop( "the left hand side of 'selection' has to contain",
         " exactly two levels (e.g. FALSE and TRUE)" )
   }
   probitDummy <- probitEndogenous == probitLevels[ 2 ]
   ## Outcome equation
   m <- match(c("formula", "data", "subset", "weights", "na.action",
                "offset"), names(mf), 0)
   mfO <- mf[c(1, m)]
   mfO$drop.unused.levels <- TRUE
   mfO$na.action <- na.pass
                                        # Here we have to keep NA-s: unobserved outcome variables may
                                        # be marked as NA-s.
   mfO[[1]] <- as.name("model.frame")
                                        # eval it as model frame
   names(mfO)[2] <- "formula"
   mfO <- eval(mfO, parent.frame())
   mtO <- attr(mfO, "terms")
   secondStepData <- model.matrix(mtO, mfO)
                                        # the explanatory variables in matrix form
   NXO <- ncol(secondStepData)
   secondStepEndogenous <- model.response(mfO, "numeric")
   # NA action
   if( !is.null( inst ) ) {
      secondStepData <- cbind( secondStepData,
         model.frame( inst, data = data, na.action = NULL ) )
   }
   firstStepOk <- rowSums( is.na( firstStepData ) ) == 0
   secondStepOk <- rowSums( is.na( secondStepData ) ) == 0
   result$dataOk <- firstStepOk & ( secondStepOk | !probitDummy )
   if( print.level > 0 ) {
      cat ( "\nEstimating 1st step Probit model . . ." )
   }
   result$probit <- probit(selection, data=data, x=TRUE, print.level=print.level - 1)
   if( print.level > 0 ) {
       cat( " OK\n" )
   }
#    data$probitLambda <- dnorm(linearPredictors(result$probit)) /
#        pnorm(linearPredictors(result$probit))
#    data$probitDelta <- data$probitLambda * ( data$probitLambda +
#                                             linearPredictors(result$probit))
   imrData <- invMillsRatio( result$probit )
   result$imrDelta <- imrData$delta1

   if( is.null( inst ) ) {
      if( print.level > 0 ) {
         cat ( "Estimating 2nd step (outcome) OLS model . . ." )
      }
      result$lm <- lm(secondStepEndogenous ~ -1 + secondStepData + imrData$IMR1,
                      subset = probitDummy )
      resid <- residuals( result$lm )
       step2coef <- coef( result$lm )
      names(step2coef) <- c(colnames(secondStepData), "invMillsRatio")
      intercept <- any(apply(model.matrix(result$lm), 2,
                             function(v) (v[1] > 0) & (all(v == v[1]))))
                                        # we have determine whether the outcome model has intercept.
                                        # This is necessary later for calculating R^2
      if( print.level > 0 ) cat( " OK\n" )
   }
   else {
      data$invMillsRatio <- imrData$IMR1
      if( print.level > 0 ) {
         cat ( "Estimating 2nd step 2SLS/IV model . . ." )
      }
      step2formula <- as.formula( paste( formula[ 2 ], "~", formula[ 3 ],
                                        "+ invMillsRatio" ) )
      formulaList <- list( step2formula )
      instImr <- as.formula( paste( "~", inst[ 2 ], "+ invMillsRatio" ) )
      library( systemfit )
      result$lm <- systemfit( "2SLS", formulaList, inst = instImr,
         data = data[ probitDummy, ] )
      resid <- residuals( result$lm )[ , 1 ]
       step2coef <- coefficients( result$lm$eq[[ 1 ]] )
      if( print.level > 0 ) cat( " OK\n" )
   }
   result$sigma <- as.numeric( sqrt( crossprod( resid ) /
      sum( probitDummy ) +
      mean(result$imrDelta[ probitDummy ] ) *
       step2coef[ "invMillsRatio" ]^2 ) )
   result$rho <-  step2coef[ "invMillsRatio" ] / result$sigma
   names(result$rho) <- NULL
                                        # otherwise the name of step2coef is left...
   result$invMillsRatio <- invMillsRatio
   if( print.level > 0 ) {
      cat ( "Calculating coefficient covariance matrix . . ." )
   }
   # the following variables are named according to Greene (2003), p. 785
   if( is.null( inst ) ) {
      xMat <- model.matrix( result$lm )
   } else {
      xMat <- result$lm$eq[[ 1 ]]$x
   }
   result$vcov <- heckitVcov( xMat,
      model.matrix( result$probit )[ probitDummy, ],
      vcov( result$probit ),
      result$rho,
      result$imrDelta[ probitDummy ],
      result$sigma )
# result$vcov <- result$sigma^2 * solve( crossprod( xMat ) ) %*%
#      ( txd2Mat %*%
#      xMat + qMat ) %*% solve( crossprod( xMat ) )
#   rm( txd2Mat, d2Vec )
   if( print.level > 0 ) cat( " OK\n" )
   result$coef <- matrix( NA, nrow = length( step2coef ), ncol = 4 )
   rownames( result$coef ) <- names( step2coef )
   colnames( result$coef ) <- c( "Estimate", "Std. Error", "t value",
      "Pr(>|t|)" )
   result$coef[ , 1 ] <- step2coef
   result$coef[ , 2 ] <- sqrt( diag( result$vcov ) )
   result$coef[ , 3 ] <- result$coef[ , 1 ] / result$coef[ , 2 ]
   result$coef[ , 4 ] <- 2 * ( 1 - pt( abs( result$coef[ , 3 ] ),
      result$lm$df ) )
   result$param <- list(index=list(betaS=seq(length=NXS), betaO=NXS + seq(length=NXO),
                        invMillsRation=NXS + NXO + 1, sigma=NXS + NXO + 2, rho=NXS + NXO + 3))
                                        # The location of results in the coef vector
   result <- c(result,
               list(outcome=list(intercept=intercept))
                                        # The 'oucome' component is intended for holding all kind of
                                        # stuff, related to the outcome equation
               )
   class( result ) <- "heckit"
   return( result )
}

summary.heckit <- function( object, ... ) {
   ## Calculate r-squared.  Note that the way lm() finds R2 is a bit naive -- it checks for intercept
   ## in the formula, but not whether the intercept is present in any of the data vectors (or matrices)
   oModel <- object$lm
   y <- model.response(model.frame(object$lm))
   if(object$outcome$intercept) {
      R2 <- sum(residuals(oModel)^2)/sum((y - mean(y))^2)
      R2adj <- 1 - (1 - R2)*(NObs(oModel) - 1)/(NObs(oModel) - NParam(oModel))
   }
   else {
      R2 <- sum(residuals(oModel)^2)/sum(y^2)
      R2adj <- 1 - (1 - R2)*(NObs(oModel))/(NObs(oModel) - NParam(oModel))
   }
   s <- c(object, r.squared=list(c(R2, R2adj)))
   class(s) <- c("summary.heckit", class(s))
   s
}

print.summary.heckit <- function( x, digits = 6, ... ) {
   Signif <- symnum( x$coef[ , 4 ], corr = FALSE, na = FALSE,
      cutpoints = c( 0, 0.001, 0.01, 0.05, 0.1, 1 ),
      symbols   = c( "***", "**", "*", "." ," " ))

   table <- cbind( round( x$coef, digits ), Signif )

   rownames( table ) <- rownames( x$coef )
   colnames( table ) <- c( "Estimate", "Std. Error", "t value",
      "Pr(>|t|)", "" )

   print( table, quote = FALSE, right = TRUE )
   cat( "---\nSignif. codes: ", attr( Signif, "legend" ), "\n" )

   if( class( x$lm ) == "lm" ) {
      rSquared <- c( summary( x$lm )$r.squared, summary( x$lm )$adj.r.squared )
   } else {
      rSquared <- c( x$lm$eq[[ 1 ]]$r2, x$lm$eq[[ 1 ]]$adjr2 )
   }
   cat( paste(
      "Multiple R-Squared:", round( rSquared[ 1 ], digits),
      "Adjusted R-Squared:", round( rSquared[ 2 ], digits),
      "\n" ) )
   cat("\n")
   invisible( x )
}
