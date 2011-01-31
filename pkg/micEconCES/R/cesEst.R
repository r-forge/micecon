cesEst <- function( yName, xNames, data, vrs = FALSE,
      method = "LM", start = NULL, lower = NULL, upper = NULL,
      rho1 = NULL, rho2 = NULL, rho = NULL, returnGridAll = FALSE, 
      random.seed = 123, rhoApprox = c( y = 5e-6, gamma = 5e-6, delta = 5e-6, 
         rho = 1e-3, nu = 5e-6 ), ... ) {

   # y = gamma * ( delta * x1^(-rho) + ( 1 - delta ) * x2^(-rho) )^(-nu/rho)
   # s = 1 / ( 1 + rho )

   checkNames( c( yName, xNames ), names( data ) )

   # abbreviated method
   if( method == "NM" ) {
      method <- "Nelder-Mead"
   }

   # number of exogenous variables
   nExog <- length( xNames )
   if( nExog == 2 ) {
      nested <- FALSE
   } else if( nExog %in% c( 3, 4 ) ) {
      nested <- TRUE
   } else {
      stop( "currently, the CES can be estimated",
         " with two inputs (normal non-nested CES)",
         " or with three or four inputs (nested CES) only" )
   }

   # obtain names of coefficients
   coefNames <- cesCoefNames( nExog = nExog, vrs = vrs, 
      returnRho1 = is.null( rho1 ), returnRho2 = is.null( rho2 ), 
      returnRho = is.null( rho ), nested = nested )

   # check rhoApprox
   if( !nested ) {
      rhoApprox <- cesCheckRhoApprox( rhoApprox = rhoApprox, withY = TRUE,
         withDeriv = TRUE )
   }

   # checking rho1
   if( !is.null( rho1 ) ) {
      if( !nested ) {
         stop( "argument 'rho1' can be used only for estimating",
            " nested CES functions" )
      } else if( !is.numeric( rho1 ) ) {
         stop( "argument 'rho1' must be either 'NULL' or numeric" )
      } else if( min( rho1 ) < -1 ) {
         stop( "the rho1s specified in argument 'rho1'",
            " must not be smaller than '-1'" )
      }
   }

   # checking rho2
   if( !is.null( rho2 ) ) {
      if( !nested ) {
         stop( "argument 'rho2' can be used only for estimating",
            " nested CES functions" )
      } else if( !is.numeric( rho2 ) ) {
         stop( "argument 'rho2' must be either 'NULL' or numeric" )
      } else if( min( rho2 ) < -1 ) {
         stop( "the rho2s specified in argument 'rho2'",
            " must not be smaller than '-1'" )
      }
   }

   # checking "rho"
   if( !is.null( rho ) ) {
      if( !is.numeric( rho ) ) {
         stop( "argument 'rho' must be either 'NULL' or numeric" )
      } else if( min( rho ) < -1 ) {
         stop( "the rhos specified in argument 'rho'",
            " must not be smaller than '-1'" )
      }
   }

   # grid search for rho_1 and rho
   if( length( rho1 ) > 1 || length( rho2 ) > 1  || length( rho ) > 1 ) {
      result <- cesEstGridRho( yName = yName, xNames = xNames,
         data = data, vrs = vrs, method = method, start = start,
         lower = lower, upper = upper,
         rho1Values = rho1, rho2Values = rho2, rhoValues = rho, 
         returnAll = returnGridAll,
         random.seed = random.seed, rhoApprox = rhoApprox, ... )
      result$call <- match.call()
      return( result )
   }
   
   # number of parameters
   nParam <- 1 + nExog + vrs + nested * 2 - ( !is.null( rho1 ) ) - 
      ( !is.null( rho2 ) ) - ( !is.null( rho ) ) - ( nested && nExog == 3 )

   # start values
   start <- cesEstStart( yName = yName, xNames = xNames, data = data,
      vrs = vrs, method = method, start = start, rho1 = rho1, rho2 = rho2, 
      rho = rho, nParam = nParam, nested = nested )

   # dertermining lower and upper bounds automatically
   if( is.null( lower ) ) {
      lower <- cesCoefBounds( vrs = vrs, returnRho1 = is.null( rho1 ), 
         returnRho2 = is.null( rho2 ), returnRho = is.null( rho ),
         method = method, lower = TRUE, nExog = nExog, nested = nested )
   }
   if( is.null( upper ) ) {
      upper <- cesCoefBounds( vrs = vrs, returnRho1 = is.null( rho1 ), 
         returnRho2 = is.null( rho2 ), returnRho = is.null( rho ),
         method = method, lower = FALSE, nExog = nExog, nested = nested )
   }

   # checking lower and upper bounds
   if( method %in% c( "L-BFGS-B", "PORT", "DE" ) ) {
      if( length( lower ) > 1 && length( lower ) != nParam ) {
         stop( "the lower bound has ", length( lower ), " elements",
            " but the model has ", nParam, " parameters" )
      }
      if( method != "DE" && any( start < lower ) ) {
         stop( "at least one starting value is smaller than its lower bound" )
      }
      if( length( upper ) > 1 && length( upper ) != nParam ) {
         stop( "the upper bound has ", length( upper ), " elements",
            " but the model has ", nParam, " parameters" )
      }
      if( method != "DE" && any( start > upper ) ) {
         stop( "at least one starting value is greater than its upper bound" )
      }
      if( length( lower ) == length( upper ) ) {
         if( any( lower > upper ) ) {
            stop( "at least one lower bound is greater than its upper bound" )
         }
      }
   } else if( max( lower ) != -Inf || min( upper ) != Inf ) {
      warning( "lower and upper bounds are ignored in method '", method, "'" )
      lower <- -Inf
      upper <- Inf
   }

   # store the (matched) call
   matchedCall <- match.call()

   # save seed of the random number generator
   if( exists( ".Random.seed" ) ) {
      savedSeed <- .Random.seed
   }

   # set seed for the random number generator (used by SANN and DE)
   set.seed( random.seed )

   # restore seed of the random number generator on exit
   # (end of function or error)
   if( exists( "savedSeed" ) ) {
      on.exit( assign( ".Random.seed", savedSeed, envir = sys.frame() ) )
   } else {
      on.exit( rm( .Random.seed, envir = sys.frame() ) )
   }

   # prepare list that will be returned
   result <- list()

   # Estimation by the Kmenta approximation
   if( method == "Kmenta" ) {
      if( nested ) {
         stop( "method 'Kmenta' is not supported for nested models" )
      }
      if( !is.null( rho ) ) {
         stop( "fixing 'rho' is currently not supported for the",
            " Kmenta approximation" )
      }
      result <- cesEstKmenta( yName = yName, xNames = xNames, data = data,
         vrs = vrs )
   } else if( method %in% c( "Nelder-Mead", "SANN", "BFGS", "CG", "L-BFGS-B" ) ) {
      if( method %in% c( "Nelder-Mead", "SANN" ) ) {
         result$optim <- optim( par = start, fn = cesRss, data = data,
            method = method, yName = yName, xNames = xNames, vrs = vrs,
            rho1 = rho1, rho2 = rho2, rho = rho, rhoApprox = rhoApprox, 
            nested = nested, ... )
      } else {
         result$optim <- optim( par = start, fn = cesRss, gr = cesRssDeriv,
            data = data, method = method, lower = lower, upper = upper, 
            yName = yName, xNames = xNames, vrs = vrs, rho1 = rho1, 
            rho2 = rho2, rho = rho, rhoApprox = rhoApprox, 
            nested = nested, ... )
      }
      result$coefficients <- result$optim$par
      result$iter <- result$optim$counts[ !is.na( result$optim$counts ) ]
      if( length( result$iter ) == 1 ) {
         result$iter <- unname( result$iter )
      }
      if( method != "SANN" ) {
         result$convergence <- result$optim$convergence == 0
      }
      result$message <- result$optim$message
      rss <- result$optim$value
   } else if( method == "LM" ) {
      # residual function
      residFun <- function( par, yName, xNames, data, vrs, rho1, rho2, rho,
            rhoApprox, nested ) {
         # add coefficients rho_1, rho_2, and rho, if they are fixed
         par <- cesCoefAddRho( coef = par, vrs = vrs, rho1 = rho1,
            rho2 = rho2, rho = rho, nExog = length( xNames ), nested = nested )
         result <- data[[ yName ]] - cesCalc( xNames = xNames,
            data = data, coef = par, rhoApprox = rhoApprox[ "y" ], 
            nested = nested )
         return( result )
      }

      # jacobian function
      jac <- function( par, yName, xNames, data, vrs, rho1, rho2, rho,
            rhoApprox, nested ) {
         # add coefficients rho_1, rho_2, and rho, if they are fixed
         par <- cesCoefAddRho( coef = par, vrs = vrs, rho1 = rho1, 
            rho2 = rho2, rho = rho, nExog = length( xNames ), nested = nested )
         return( -c( cesDerivCoef( par = par, xNames = xNames, data = data,
            vrs = vrs, returnRho1 = is.null( rho1 ), 
            returnRho2 = is.null( rho2 ), returnRho = is.null( rho ),
            rhoApprox = rhoApprox, nested = nested ) ) )
      }

      # perform fit
      result$nls.lm <- nls.lm( par = start, fn = residFun, data = data,
         jac = jac, yName = yName, xNames = xNames, vrs = vrs, rho1 = rho1, 
         rho2 = rho2, rho = rho, rhoApprox = rhoApprox, nested = nested, ... )
      result$coefficients <- result$nls.lm$par
      result$iter <- result$nls.lm$niter
      result$convergence <- result$nls.lm$info > 0 && result$nls.lm$info < 5
      result$message <- result$nls.lm$message
      rss <- result$nls.lm$deviance
   } else if( method == "Newton" ) {
      cesRss2 <- function( par, yName, xNames, data, vrs, rho1, rho2, rho, 
            rhoApprox, nested ) {
         result <- cesRss( par = par, yName = yName, xNames = xNames,
            data = data, vrs = vrs, rho1 = rho1, rho2 = rho2, rho = rho, 
            rhoApprox = rhoApprox[ "y" ], nested = nested )
         attributes( result )$gradient <- cesRssDeriv( par = par, 
            yName = yName, xNames = xNames, data = data, vrs = vrs, rho1 = rho1,
            rho2 = rho2, rho = rho, rhoApprox = rhoApprox, nested = nested )
         return( result )
      }
      # save current setting for warning messages and suppress warning messages
      warnSaved <- options()$warn
      options( warn = -1 )
      # perform fit
      result$nlm <- nlm( f = cesRss2, p = start, data = data,
         yName = yName, xNames = xNames, vrs = vrs, rho1 = rho1, rho2 = rho2, 
         rho = rho, rhoApprox = rhoApprox, nested = nested, ... )
      # restore previous setting for warning messages
      options( warn = warnSaved )
      # extract results
      result$coefficients <- result$nlm$estimate
      result$iter <- result$nlm$iterations
      result$convergence <- result$nlm$code <= 2
      rss <- result$nlm$minimum
   } else if( method == "PORT" ) {
      result$nlminb <- nlminb( start = start, objective = cesRss,
         gradient = cesRssDeriv, data = data, yName = yName, xNames = xNames,
         vrs = vrs, rho1 = rho1, rho2 = rho2, rho = rho, 
         lower = lower, upper = upper, 
         rhoApprox = rhoApprox, nested = nested, ... )
      result$coefficients <- result$nlminb$par
      result$iter <- result$nlminb$iterations
      result$convergence <- result$nlminb$convergence == 0
      result$message <- result$nlminb$message
      rss <- result$nlminb$objective
   } else if( method == "DE" ) {
      result$DEoptim <- DEoptim( fn = cesRss, lower = lower,
         upper = upper, data = data, yName = yName, xNames = xNames,
         vrs = vrs, rho1 = rho1, rho2 = rho2, rho = rho, 
         rhoApprox = rhoApprox, nested = nested, ... )
      result$coefficients <- result$DEoptim$optim$bestmem
      result$iter <- result$DEoptim$optim$iter
      rss <- result$DEoptim$optim$bestval
   } else if( method == "nls" ) {
      if( !is.null( rho1 ) ) {
         warning( "ignoring argument 'rho1'" )
      }
      if( !is.null( rho2 ) ) {
         warning( "ignoring argument 'rho2'" )
      }
      if( !is.null( rho ) ) {
         warning( "ignoring argument 'rho'" )
      }
      if( nested && nExog == 3 ) {
         nlsFormula <- as.formula( paste( yName,
            " ~ gamma * ( delta_2 * ",
            "( delta_1 * ", xNames[ 1 ], "^(-rho_1)",
            " + ( 1 - delta_1 ) * ", xNames[ 2 ], "^(-rho_1) )",
            "^( rho / rho_1 ) +",
            " ( 1 - delta_2 ) * ", xNames[ 3 ], "^(-rho) )",
            "^( - ", ifelse( vrs, "nu", "1" ), " / rho )",
            sep = "" ) )
      } else if( nested && nExog == 4 ) {
         nlsFormula <- as.formula( paste( yName,
            " ~ gamma * ( delta_3 * ( delta_1 * ", xNames[ 1 ], "^(-rho_1)",
            " + ( 1 - delta_1 ) * ", xNames[ 2 ], "^(-rho_1) )",
            "^( rho / rho_1 ) +",
            " ( 1 - delta_3 ) * ( delta_2 * ", xNames[ 3 ], "^(-rho_2)",
            " + ( 1 - delta_2 ) * ", xNames[ 4 ], "^(-rho_2) )",
            "^( rho / rho_2 ) )",
            "^( - ", ifelse( vrs, "nu", "1" ), " / rho )",
            sep = "" ) )
      } else {
         nlsFormula <- as.formula( paste( yName,
            " ~ gamma * ( delta * ", xNames[ 1 ], "^(-rho)",
            " + ( 1 - delta ) * ", xNames[ 2 ], "^(-rho) )",
            "^( - ", ifelse( vrs, "nu", "1" ), " / rho )",
            sep = "" ) )
      }
      result$nls <- nls( formula = nlsFormula, data = data, start = start,
         ... )
      result$coefficients <- coef( result$nls )
      result$iter <- result$nls$convInfo$finIter
      result$convergence <- result$nls$convInfo$isConv
      if( result$nls$convInfo$stopMessage != "converged" ) {
         result$message <- result$nls$convInfo$stopMessage
      }
      rss <- deviance( result$nls )
   } else {
      stop( "argument 'method' must be either 'Nelder-Mead', 'BFGS',",
         " 'CG', 'L-BFGS-B', 'SANN', 'LM', 'Newton', 'PORT',",
         " 'DE', 'nls', or 'Kmenta'" )
   }

   # add names to estimated coefficients
   names( result$coefficients ) <- coefNames

   # add rho_1, rho_2, and rho, if they are fixed
   result$coefficients <- cesCoefAddRho( coef = result$coefficients,
      vrs = vrs, rho1 = rho1, rho2 = rho2, rho = rho, 
      nExog = nExog, nested = nested )


   # return also the call
   result$call <- matchedCall

   # return the method used for the estimation
   result$method <- method

   # return the starting values
   result$start <- start

   # return lower and upper bounds
   result$lower <- lower
   result$upper <- upper

   # return fixed rho_1
   result$rho1 <- rho1

   # return fixed rho_2
   result$rho2 <- rho2

   # return fixed rho
   result$rho <- rho

   # fitted values
   result$fitted.values <- cesCalc( xNames = xNames, data = data,
      coef = result$coefficients, rhoApprox = rhoApprox[ "y" ], nested = nested )

   # residuals
   result$residuals <- data[[ yName ]] - result$fitted.values

   # sum of squared residuals
   result$rss <- sum( result$residuals^2 )
   if( method != "Kmenta" ) {
      if( !isTRUE( all.equal( rss, result$rss, tolerance = 1e-3 ) ) ) {
         warning( "internal problem: the minimum of the objective function",
            " returned by the solver (", rss, ") is not equal to the",
            " RSS calculated from the residuals (", result$rss, ")" )
      }
   }

   # unscaled covariance matrix
   gradients <- cesDerivCoef( par = result$coefficients, xNames = xNames,
      data = data, vrs = vrs, rhoApprox = rhoApprox, nested = nested )
   result$cov.unscaled <- try( chol2inv( chol( crossprod( gradients ) ) ),
      silent = TRUE )
   if( !is.matrix( result$cov.unscaled ) ) {
      result$cov.unscaled <- matrix( NA, nrow = length( result$coefficients ),
         ncol = length( result$coefficients ) )
   }
   rownames( result$cov.unscaled ) <- names( result$coefficients )
   colnames( result$cov.unscaled ) <- names( result$coefficients )

   class( result ) <- c( "cesEst", class( result ) )
   return( result )
}

