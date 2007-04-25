## heckit
vcov.selection <- function( object, part = "full", ... ) {
   if( object$method == "ml" ){
      result <- vcov.maxlik( object, ... )
   } else if( object$method == "2step" ) {
      if( object$tobitType == 2 ) {
         result <- vcov.heckit( object, part = part, ... )
      } else if( object$tobitType == 5 ) {
         result <- vcov.heckit5( object, part = part, ... )
      }
   }
   return( result )
}

## heckit
vcov.heckit <- function(object, part="outcome", ...) {
   if(part=="outcome") {
      i <- c(object$param$index$betaO, object$param$index$Mills)
      vc <- object$vcov[i,i]
   }
   else if(part=="full") {
      vc <- object$vcov
   }
   else
       stop("'part' must be either 'outcome' or 'full'")
   vc
  }

## heckit 5
vcov.heckit5 <- function(object, part="outcome", ...) {
   if(part=="outcome") {
      i <- c(object$param$index$betaO1, object$param$index$invMillsRatio1,
             object$param$index$betaO2, object$param$index$invMillsRatio2)
      vc <- object$vcov[i,i]
   }
   else if(part=="full") {
      vc <- object$vcov
   }
   else
       stop("'part' must be either 'outcome' or 'full'")
   vc
}
