durbinWatsonTest.cesEst <- function( model, ... ) {
   if( is.null( model$grad ) ) {
      stop( "the model object must have a component 'grad';",
         " please re-estimate the model with argument 'returnGrad'",
         " set to 'TRUE'" )
   }
   gMat <- model$grad
   mCoef <- coef( model )
   yDW <- residuals( model ) + gMat %*% mCoef
   dwReg <- lm( yDW ~ gMat - 1 )
   result <- durbinWatsonTest( dwReg, ... )
   return( result )
}