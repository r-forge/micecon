cesEst <- function( yName, xNames, data, vrs = FALSE,
      method = "Nelder-Mead", startVal = NULL, lower = NULL, upper = NULL,
      ... ) {

   # y = gamma * ( delta * x1^(-rho) + ( 1 - delta ) * x2^(-rho) )^(-phi/rho)
   # s = 1 / ( 1 + rho )

   checkNames( c( yName, xNames ), names( data ) )

   if( length( xNames ) != 2 ) {
      stop( "currently, argument 'xNames' must contain exactly",
         " two variable names" )
   }

   # number of parameters
   nParam <- 3 + vrs

   # start values
   if( method %in% c( "Kmenta", "DE" ) ) {
      if( !is.null( startVal ) ) {
         warning( "ignoring starting values because they are not required",
            " for method '", method, "'" )
      }
   } else {
      if( is.null( startVal ) ) {
         startVal <- c( 1, 0.5, 0.25, 1 )[ 1:( 3 + vrs ) ]
         yTemp <- cesCalc( xNames = xNames, data = data, coef = startVal )
         startVal[ 1 ] <- mean( data[[ yName ]], na.rm = TRUE ) /
            mean( yTemp, na.rm = TRUE )
      }
      if( length( startVal ) != 4 && vrs ) {
         stop( "a CES function with 2 explanatory variables and VRS has 4",
            " parameters but you provided ", length( startVal ),
            " starting values" )
      } else if( length( startVal ) != 3 && !vrs ) {
         stop( "a CES function with 2 explanatory variables and CRS has 3",
            " parameters but you provided ", length( startVal ),
            " starting values" )
      }
      names( startVal ) <- c( "gamma", "delta", "rho", "phi" )[ 1:( 3 + vrs ) ]
   }

   # dertermining lower and upper bounds automatically
   if( is.null( lower ) ) {
      if( method %in% c( "L-BFGS-B", "PORT", "DE" ) ) {
         lower <- c( 0, 0, -1, 0 )[ 1:nParam ]
      } else {
         lower <- -Inf
      }
   }
   if( is.null( upper ) ) {
      if( method %in% c( "L-BFGS-B", "PORT" ) ) {
         upper <- c( Inf, 1, Inf, Inf )[ 1:nParam ]
      } else if( method == "DE" ) {
         upper <- c( 1e10, 1, 10, 10 )[ 1:nParam ]
      } else {
         upper <- Inf
      }
   }

   # checking lower and upper bounds
   if( method %in% c( "L-BFGS-B", "PORT", "DE" ) ) {
      if( length( lower ) > 1 && length( lower ) != nParam ) {
         stop( "the lower bound has ", length( lower ), " elements",
            " but the model has ", nParam, " parameters" )
      }
      if( method != "DE" && any( startVal < lower ) ) {
         stop( "at least one starting value is smaller than its lower bound" )
      }
      if( length( upper ) > 1 && length( upper ) != nParam ) {
         stop( "the upper bound has ", length( upper ), " elements",
            " but the model has ", nParam, " parameters" )
      }
      if( method != "DE" && any( startVal > upper ) ) {
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

   # prepare data for estimation
   estData <- data.frame( y = data[[ yName ]],
      x1 = data[[ xNames[ 1 ] ]], x2 = data[[ xNames[ 2 ] ]] )

   # prepare list that will be returned
   result <- list()

   # Estimation by the Kmenta approximation
   if( method == "Kmenta" ) {
      result <- cesEstKmenta( yName = yName, xNames = xNames, data = data,
         vrs = vrs )
   } else if( method %in% c( "Nelder-Mead", "SANN", "BFGS", "CG", "L-BFGS-B" ) ) {
      if( method %in% c( "Nelder-Mead", "SANN" ) ) {
         result$optim <- optim( par = startVal, fn = cesRss, data = estData,
            method = method, ... )
      } else {
         result$optim <- optim( par = startVal, fn = cesRss, gr = cesRssDeriv,
            data = estData, method = method, lower = lower, upper = upper, ... )
      }
      result$coefficients <- result$optim$par
      result$iter <- result$optim$counts[ !is.na( result$optim$counts ) ]
      if( length( result$iter ) == 1 ) {
         result$iter <- unname( result$iter )
      }
      result$convergence <- result$optim$convergence == 0
      result$message <- result$optim$message
   } else if( method == "LM" ) {
      # residual function
      residFun <- function( par, data2 ) {
         result <- data2$y - cesCalc( xNames = c( "x1", "x2" ),
            data = data2, coef = par )
         return( result )
      }

      # jacobian function
      jac <- function( par, data3 ) {
         return( -c( cesDerivCoef( par = par, data = data3 ) ) )
      }

      # perform fit
      result$nls.lm <- nls.lm( par = startVal, fn = residFun, data = estData,
         jac = jac, ... )
      result$coefficients <- result$nls.lm$par
      result$iter <- result$nls.lm$niter
      result$convergence <- result$nls.lm$info > 0 && result$nls.lm$info < 5
      result$message <- result$nls.lm$message
   } else if( method == "Newton" ) {
      cesRss2 <- function( par, data ) {
         result <- cesRss( par = par, data = data )
         attributes( result )$gradient <- cesRssDeriv( par = par, data = data )
         return( result )
      }
      # save current setting for warning messages and suppress warning messages
      warnSaved <- options()$warn
      options( warn = -1 )
      # perform fit
      result$nlm <- nlm( f = cesRss2, p = startVal, data = estData, ... )
      # restore previous setting for warning messages
      options( warn = warnSaved )
      # extract results
      result$coefficients <- result$nlm$estimate
      result$iter <- result$nlm$iterations
      result$convergence <- result$nlm$code <= 2
   } else if( method == "PORT" ) {
      result$nlminb <- nlminb( start = startVal, objective = cesRss,
         gradient = cesRssDeriv, data = estData,
         lower = lower, upper = upper, ... )
      result$coefficients <- result$nlminb$par
      result$iter <- result$nlminb$iterations
      result$convergence <- result$nlminb$convergence == 0
      result$message <- result$nlminb$message
   } else if( method == "DE" ) {
      result$DEoptim <- DEoptim( fn = cesRss, lower = lower,
         upper = upper, data = estData, ... )
      result$coefficients <- result$DEoptim$optim$bestmem
      names( result$coefficients ) <-
         c( "gamma", "delta", "rho", "phi" )[ 1:( 3 + vrs ) ]
      result$iter <- result$DEoptim$optim$iter
   } else {
      stop( "argument 'method' must be either 'Nelder-Mead', 'BFGS',",
         " 'CG', 'L-BFGS-B', 'SANN', 'LM', 'Newton', 'PORT',",
         " 'DE', or 'Kmenta'" )
   }

   # return also the call
   result$call <- matchedCall

   # return the method used for the estimation
   result$method <- method

   # return the starting values
   if( ! method %in% c( "Kmenta", "DE" ) ) {
      result$startVal <- startVal
   }

   # return lower and upper bounds
   result$lower <- lower
   result$upper <- upper

   # fitted values
   result$fitted.values <- cesCalc( xNames = xNames, data = data,
      coef = result$coefficients )

   # residuals
   result$residuals <- estData$y - result$fitted.values

   # unscaled covariance matrix
   gradients <- cesDerivCoef( par = result$coefficients, data = estData )
   result$cov.unscaled <- try( chol2inv( chol( crossprod( gradients ) ) ),
      silent = TRUE )
   if( !is.matrix( result$cov.unscaled ) ) {
      result$cov.unscaled <- matrix( NA, nrow = length( result$coefficients ),
         ncol = length( result$coefficients ) )
   }
   rownames( result$cov.unscaled ) <- names( result$coefficients )
   colnames( result$cov.unscaled ) <- names( result$coefficients )

   # nonlinear least squares
#    result$nls <- nls(
#       y ~ gamma * ( delta * x1^rho + ( 1 - delta ) * x2^rho )^(1/rho),
#       data = estData, start = result$startVal, trace = TRUE,
#       algorithm = "port", lower = c( -Inf, 0.01 , -Inf ),
#       upper = c( Inf, 0.99, Inf ) )

   class( result ) <- c( "cesEst", class( result ) )
   return( result )
}

