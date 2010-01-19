cesEst <- function( yName, xNames, data, vrs = FALSE,
      method = "Nelder-Mead", startVal = NULL, lower = NULL, upper = NULL,
      rho = NULL, ... ) {

   # y = gamma * ( delta * x1^(-rho) + ( 1 - delta ) * x2^(-rho) )^(-nu/rho)
   # s = 1 / ( 1 + rho )

   checkNames( c( yName, xNames ), names( data ) )

   # number of exogenous variables
   nExog <- length( xNames )
   if( nExog != 2 ) {
      stop( "currently, argument 'xNames' must contain exactly",
         " two variable names" )
   }

   # checking "rho"
   if( !is.null( rho ) ) {
      if( !is.numeric( rho ) || length( rho ) != 1 ) {
         stop( "argument 'rho' must be either 'NULL' or a numeric scalar" )
      } else if( rho < -1 ) {
         stop( "argument 'rho' must not be smaller than '-1'" )
      }
   }

   # number of parameters
   nParam <- 3 + vrs - length( rho )

   # start values
   if( method %in% c( "Kmenta", "DE" ) ) {
      if( !is.null( startVal ) ) {
         warning( "ignoring starting values because they are not required",
            " for method '", method, "'" )
         startVal <- NULL
      }
   } else {
      if( is.null( startVal ) ) {
         rhoStart <- ifelse( is.null( rho ), 0.25, rho )
         startVal <- c( 1, 0.5, rhoStart, 1 )[ 1:( 3 + vrs ) ]
         yTemp <- cesCalc( xNames = xNames, data = data, coef = startVal )
         startVal[ 1 ] <- mean( data[[ yName ]], na.rm = TRUE ) /
            mean( yTemp, na.rm = TRUE )
         if( !is.null( rho ) ) {
            startVal <- startVal[ -3 ]
         }
      }
      if( length( startVal ) != nParam ) {
         stop( "wrong number of starting values:",
            " you provided ", length( startVal ), " values",
            " but the model has ", nParam, " parameters" )
      }
      names( startVal ) <- cesCoefNames( nExog, vrs, returnRho = is.null( rho ) )
   }

   # dertermining lower and upper bounds automatically
   if( is.null( lower ) ) {
      lower <- cesCoefBounds( vrs = vrs, returnRho = is.null( rho ),
         method = method, lower = TRUE )
   }
   if( is.null( upper ) ) {
      upper <- cesCoefBounds( vrs = vrs, returnRho = is.null( rho ),
         method = method, lower = FALSE )
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

   # prepare list that will be returned
   result <- list()

   # Estimation by the Kmenta approximation
   if( method == "Kmenta" ) {
      if( !is.null( rho ) ) {
         stop( "fixing 'rho' is currently not supported for the",
            " Kmenta approximation" )
      }
      result <- cesEstKmenta( yName = yName, xNames = xNames, data = data,
         vrs = vrs )
   } else if( method %in% c( "Nelder-Mead", "SANN", "BFGS", "CG", "L-BFGS-B" ) ) {
      if( method %in% c( "Nelder-Mead", "SANN" ) ) {
         result$optim <- optim( par = startVal, fn = cesRss, data = data,
            method = method, yName = yName, xNames = xNames, vrs = vrs,
            rho = rho, ... )
      } else {
         result$optim <- optim( par = startVal, fn = cesRss, gr = cesRssDeriv,
            data = data, method = method, lower = lower, upper = upper, 
            yName = yName, xNames = xNames, vrs = vrs, rho = rho, ... )
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
      residFun <- function( par, yName, xNames, data, vrs, rho ) {
         # add coefficient 'rho' if it is fixed
         par <- cesCoefAddRho( coef = par, vrs = vrs, rho = rho )
         result <- data[[ yName ]] - cesCalc( xNames = xNames,
            data = data, coef = par )
         return( result )
      }

      # jacobian function
      jac <- function( par, yName, xNames, data, vrs, rho ) {
         # add coefficient 'rho' if it is fixed
         par <- cesCoefAddRho( coef = par, vrs = vrs, rho = rho )
         return( -c( cesDerivCoef( par = par, xNames = xNames, data = data,
            vrs = vrs, returnRho = is.null( rho ) ) ) )
      }

      # perform fit
      result$nls.lm <- nls.lm( par = startVal, fn = residFun, data = data,
         jac = jac, yName = yName, xNames = xNames, vrs = vrs, rho = rho, ... )
      result$coefficients <- result$nls.lm$par
      result$iter <- result$nls.lm$niter
      result$convergence <- result$nls.lm$info > 0 && result$nls.lm$info < 5
      result$message <- result$nls.lm$message
   } else if( method == "Newton" ) {
      cesRss2 <- function( par, yName, xNames, data, vrs, rho ) {
         result <- cesRss( par = par, yName = yName, xNames = xNames,
            data = data, vrs = vrs, rho = rho )
         attributes( result )$gradient <- cesRssDeriv( par = par, 
            yName = yName, xNames = xNames, data = data, vrs = vrs, rho = rho )
         return( result )
      }
      # save current setting for warning messages and suppress warning messages
      warnSaved <- options()$warn
      options( warn = -1 )
      # perform fit
      result$nlm <- nlm( f = cesRss2, p = startVal, data = data,
         yName = yName, xNames = xNames, vrs = vrs, rho = rho, ... )
      # restore previous setting for warning messages
      options( warn = warnSaved )
      # extract results
      result$coefficients <- result$nlm$estimate
      result$iter <- result$nlm$iterations
      result$convergence <- result$nlm$code <= 2
   } else if( method == "PORT" ) {
      result$nlminb <- nlminb( start = startVal, objective = cesRss,
         gradient = cesRssDeriv, data = data, yName = yName, xNames = xNames,
         vrs = vrs, rho = rho, lower = lower, upper = upper, ... )
      result$coefficients <- result$nlminb$par
      result$iter <- result$nlminb$iterations
      result$convergence <- result$nlminb$convergence == 0
      result$message <- result$nlminb$message
   } else if( method == "DE" ) {
      result$DEoptim <- DEoptim( fn = cesRss, lower = lower,
         upper = upper, data = data, yName = yName, xNames = xNames,
         vrs = vrs, rho = rho, ... )
      result$coefficients <- result$DEoptim$optim$bestmem
      names( result$coefficients ) <- cesCoefNames( nExog, vrs,
         returnRho = is.null( rho ) )
      result$iter <- result$DEoptim$optim$iter
   } else {
      stop( "argument 'method' must be either 'Nelder-Mead', 'BFGS',",
         " 'CG', 'L-BFGS-B', 'SANN', 'LM', 'Newton', 'PORT',",
         " 'DE', or 'Kmenta'" )
   }

   # add the 'rho' if it is fixed
   result$coefficients <- cesCoefAddRho( coef = result$coefficients,
      vrs = vrs, rho = rho )

   # return also the call
   result$call <- matchedCall

   # return the method used for the estimation
   result$method <- method

   # return the starting values
   result$startVal <- startVal

   # return lower and upper bounds
   result$lower <- lower
   result$upper <- upper

   # return fixed rho
   result$rho <- rho

   # fitted values
   result$fitted.values <- cesCalc( xNames = xNames, data = data,
      coef = result$coefficients )

   # residuals
   result$residuals <- data[[ yName ]] - result$fitted.values

   # sum of squared residuals
   result$rss <- sum( result$residuals^2 )

   # unscaled covariance matrix
   gradients <- cesDerivCoef( par = result$coefficients, xNames = xNames,
      data = data, vrs = vrs )
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

