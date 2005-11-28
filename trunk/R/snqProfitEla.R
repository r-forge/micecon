## ===== calculation of elasticities from beta matrix ===
snqProfitEla <- function( beta, prices, quant, weights, 
   coefVcov = NULL, df = NULL ) {
   if( !is.matrix( beta ) ) {
      stop( "argument 'beta' must be a matrix" )
   }
   if( nrow( beta ) != ncol( beta ) ) {
      stop( "argument 'beta' must be a quadratic matrix" )
   }
   if( length( prices ) != length( quant ) ) {
      stop( "arguments 'prices' and 'quant' must have the same length" )
   }
   if( length( prices ) != length( weights ) ) {
      stop( "arguments 'prices' and 'weights' must have the same length" )
   }
   if( nrow( beta ) != length( prices ) ) {
      stop( "arguments 'prices' must have as many elements as",
         " argument 'beta' has rows" )
   }
   nNetput  <- ncol( beta )
   prices   <- unlist( prices )
   quant    <- unlist( quant )
   hessian  <- snqProfitHessian( beta, prices, weights )
   result   <- list()
   result$ela <- hessian * array( 1, c( nNetput ) ) %*% t( prices ) /
                  quant %*% t( array( 1, c( nNetput ) ) )
   if( !is.null( names( quant ) ) ) {
      rownames( result$ela ) <- names( quant )
   } else {
      rownames( result$ela ) <- paste( "q", 1:nNetput, sep = "" )
   }
   if( !is.null( names( prices ) ) ) {
      colnames( result$ela ) <- names( prices )
   } else {
      colnames( result$ela ) <- paste( "p", 1:nNetput, sep = "" )
   }
   if( !is.null( coefVcov ) ) {
      jacobian <- snqProfitElaJacobian( beta, prices, quant, weights ) 
      betaIndex <- grep( "beta", rownames( coefVcov ) )
      betaVcov <- coefVcov[ betaIndex, betaIndex ]
      result$vcov <- jacobian %*% betaVcov %*% t( jacobian )
      result$stEr <- matrix( diag( result$vcov )^0.5, nrow = nNetput,
         byrow = TRUE )
      result$tval <- result$ela / result$stEr
      if( !is.null( df ) ) {
         result$pval <- 2 * pt( abs( result$tval ), df, 
            lower.tail = FALSE )
      }
   }
   class( result ) <- "snqProfitEla"
   return( result )
}
