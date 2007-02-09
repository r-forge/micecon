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
   m <- match(c("selection", "data", "subset", "weights", 
                "offset"), names(mf), 0)
   mfS <- mf[c(1, m)]
   mfS$na.action <- na.pass
   mfS$drop.unused.levels <- TRUE
   mfS[[1]] <- as.name("model.frame")
   names(mfS)[2] <- "formula"
                                        # model.frame requires the parameter to
                                        # be 'formula'
   mfS <- eval(mfS, parent.frame())
   badRowS <- apply(mfS, 1, function(v) any(is.na(v)))
                                        # check for NA-s.  Because we have to find NA-s in several
                                        # frames, we cannot use the standard na.* functions here.
                                        # Find bad rows and remove them later.
   mtS <- attr(mfS, "terms")
   XS <- model.matrix(mtS, mfS)
   NXS <- ncol(XS)
   YS <- factor(model.response(mfS, "numeric"))
   probitLevels <- levels( as.factor( YS ) )
   if( length( probitLevels ) != 2 ) {
      stop( "the left hand side of 'selection' has to contain",
         " exactly two levels (e.g. FALSE and TRUE)" )
   }
   probitDummy <- YS == probitLevels[ 2 ]
   ## Outcome equation
   m <- match(c("formula", "data", "subset", "weights",
                "offset"), names(mf), 0)
   mfO <- mf[c(1, m)]
   mfO$na.action <- na.pass
   mfO$drop.unused.levels <- TRUE
   mfO$na.action <- na.pass
                                        # Here we have to keep NA-s: unobserved outcome variables may
                                        # be marked as NA-s.
   mfO[[1]] <- as.name("model.frame")
                                        # eval it as model frame
   names(mfO)[2] <- "formula"
   mfO <- eval(mfO, parent.frame())
   mtO <- attr(mfO, "terms")
   XO <- model.matrix(mtO, mfO)
                                        # the explanatory variables in matrix form
   NXO <- ncol(XO)
   YO <- model.response(mfO, "numeric")
   ## Remove NA observations
   badRowO <- apply(mfO, 1, function(v) any(is.na(v))) & (!is.na(YS) &YS==1)
   badRow <- badRowS | badRowO
                                        # rows in outcome, which contain NA and are observable -> bad too
   if(print.level > 0) {
      cat(sum(badRow), "invalid observations\n")
   }
   XS <- XS[!badRow,]
   YS <- YS[!badRow]
   XO <- XO[!badRow,]
   YO <- YO[!badRow]
   ##
   NParam <- NXS + NXO + 3
                                        # invMillsRation, sigma, rho = 3
   if( print.level > 0 ) {
      cat ( "\nEstimating 1st step Probit model . . ." )
   }
   result$probit <- probit(selection, data=data, x=TRUE, print.level=print.level - 1,
                           subset=eval(!badRow))
   if( print.level > 0 ) {
       cat( " OK\n" )
   }
   imrData <- invMillsRatio( result$probit )
   result$imrDelta <- imrData$delta1
   ## ---- Outcome estimation -----
   if( is.null( inst ) ) {
      if( print.level > 0 ) {
         cat ( "Estimating 2nd step (outcome) OLS model . . ." )
      }
      result$lm <- lm(YO ~ -1 + XO + imrData$IMR1,
                      subset = probitDummy )
      intercept <- any(apply(model.matrix(result$lm), 2,
                             function(v) (v[1] > 0) & (all(v == v[1]))))
                                        # we have determine whether the outcome model has intercept.
                                        # This is necessary later for calculating R^2
      resid <- residuals( result$lm )
      step2coef <- coef( result$lm )
      names(step2coef) <- c(colnames(XO), "invMillsRatio")
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
      intercept = FALSE
                                        # we calculate R^2 differently here (hopefully)
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
   coefs <- c(coef(result$probit), step2coef, sigma=result$sigma, rho=result$rho)
   if( print.level > 0 ) {
      cat ( "Calculating coefficient covariance matrix . . ." )
   }
                                        # the following variables are named according to Greene (2003), p. 785
   if( is.null( inst ) ) {
      xMat <- model.matrix( result$lm )
   } else {
      xMat <- result$lm$eq[[ 1 ]]$x
   }
   ## Now indices for packing the separate outcomes into full outcome vectors
   iBetaS <- seq(length=NXS)
   iBetaO <- NXS + seq(length=NXO)
   iMills <- NXS + NXO + 1
   iSigma <- iMills + 1
   iRho <- iSigma + 1
   ## Varcovar matrix.  Fill only a few parts, rest will remain NA
   vc <- matrix(0, NParam, NParam)
   colnames(vc) <- row.names(vc) <- names(coefs)
   vc[] <- NA
   vc[iBetaS,iBetaS] <- vcov(result$probit)
   vc[c(iBetaO,iMills), c(iBetaO,iMills)] <- heckitVcov( xMat,
                                                        model.matrix( result$probit )[ probitDummy, ],
                                                        vcov( result$probit ),
                                                        result$rho,
                                                        result$imrDelta[ probitDummy ],
                                                        result$sigma )
   result$vcov <- vc
   ##
   if( print.level > 0 )
       cat( " OK\n" )
   result$coefficients <- coefs
                                        # for coef() etc. methods
   ## the 'param' component is intended to all kind of technical info
   result$param <- list(index=list(betaS=iBetaS, betaO=iBetaO, invMillsRatio=iMills,
                                   sigma=iSigma, rho=iRho),
                                        # The location of results in the coef vector
                        oIntercept=intercept,
                        NParam=NParam)
   class( result ) <- "heckit"
   return( result )
}

summary.heckit <- function( object, part="outcome", ... ) {
   ## Calculate r-squared.  Note that the way lm() finds R2 is a bit naive -- it checks for intercept
   ## in the formula, but not whether the intercept is present in any of the data vectors (or matrices)
   oModel <- object$lm
   if(class(oModel) == "lm") {
      y <- model.response(model.frame(object$lm))
      if(object$param$oIntercept) {
         R2 <- 1 - sum(residuals(oModel)^2)/sum((y - mean(y))^2)
         R2adj <- 1 - (1 - R2)*(NObs(oModel) - 1)/(NObs(oModel) - NParam(oModel))
      }
      else {
         R2 <- 1 - sum(residuals(oModel)^2)/sum(y^2)
         R2adj <- 1 - (1 - R2)*(NObs(oModel))/(NObs(oModel) - NParam(oModel))
      }
   }
   else {
      R2 <- object$lm$eq[[ 1 ]]$r2
      R2adj <- object$lm$eq[[ 1 ]]$adjr2
   }
   if(part=="full") {
      i <- c(object$param$index$betaS, object$param$index$betaO, object$param$invMillsRatio)
   }
   else if(part=="outcome") {
      i <- object$param$index$betaO
   }
   coefs <- coef(object, part="full")[i]
   coefficients <- matrix(0, nrow = length(coefs), ncol = 4 )
   rownames(coefficients) <- names(coefs)
   colnames(coefficients) <- c( "Estimate", "Std. Error", "t value", "Pr(>|t|)" )
   coefficients[ , 1 ] <- coefs
   coefficients[ , 2 ] <- sqrt( diag(vcov(object, part="full")[i,i]))
   coefficients[ , 3 ] <- coefs/coefficients[, 2 ]
   coefficients[ , 4 ] <- 2*pt(abs(coefficients[,3]), object$lm$df, lower.tail=FALSE)
   object$coefficients <- coefficients
   s <- c(object,
          rSquared=list(c(R2, R2adj)))
   class(s) <- c("summary.heckit", class(s))
   s
}

print.summary.heckit <- function( x,
                                 digits=max(3, getOption("digits") - 3),
                                 signif.stars=getOption("show.signif.stars"),
                                 ...) {
   cat("Estimates:\n")
   printCoefmat(x$coefficients)
   cat("Multiple R-Squared:", round(x$rSquared[ 1 ], digits),
       ",\tAdjusted R-Squared:", round(x$rSquared[ 2 ], digits), "\n", sep="")
   cat("Residual correlation: ", x$rho, ", variance: ", x$sigma, "\n", sep="")
   invisible( x )
}
