translogCheckCurvature <- function( xNames, data, coef, convexity = TRUE,
   quasi = FALSE, quadHalf = TRUE, dataLogged = FALSE, ... ) {

   result <- list()

   hessian <- translogHessian( xNames = xNames, data = data, coef = coef,
      quadHalf = quadHalf, dataLogged = dataLogged, bordered = quasi )

   if( quasi ) {
      if( convexity ) {
         result$obs <- quasiconvexity( hessian, ... )
      } else {
         result$obs <- quasiconcavity( hessian, ... )
      }
   } else {
      semidef <- semidefiniteness( hessian, ... )
      if( convexity ) {
         result$obs <- semidef$positive
      } else {
         result$obs <- semidef$negative
      }
   }

   result$convexity <- convexity
   result$quasi     <- quasi

   class( result ) <- "translogCheckCurvature"
   return( result )
}
