## probit
vcov.probit <- function(object, ...) {
  result <- vcov.maxLik( object )
  if(!is.null(result))
      rownames( result ) <- colnames( result ) <- names( object$estimate )
  return( result )
}

## heckit
vcov.selection <- function( object, ... ) {
   if( object$method == "ml" ){
      result <- vcov.maxlik( object, ... )
   } else if( object$method == "2step" ) {
      if( object$tobitType == 2 ) {
         result <- vcov.heckit( object, ... )
      } else if( object$tobitType == 5 ) {
         result <- vcov.heckit5( object, ... )
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

## maxLik
vcov.maxLik <- function(object, ...) {
   ## if exists $varcovar, take it
   if(!is.null(object$varcovar))
       return(object$varcovar)
   ## otherwise invert hessian
   activePar <- activePar(object)
   if(min(abs(eigen(Hessian(object)[activePar,activePar],
                    symmetric=TRUE, only.values=TRUE)$values)) > 1e-6) {
      varcovar <- matrix(0, NParam(object), NParam(object))
      varcovar[activePar,activePar] <- solve(-Hessian(object)[activePar,activePar])
   }
   else
       varcovar <- NULL
   varcovar
}
