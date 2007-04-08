## probit
vcov.probit <- function(object, ...) {
  result <- vcov.maxLik( object )
  if(!is.null(result))
      rownames( result ) <- colnames( result ) <- names( object$estimate )
  return( result )
}

## heckit
vcov.heckit <- function(object, part="outcome", ...) {
   if(part=="outcome") {
      i <- c(object$param$index$betaO, object$param$index$invMillsRatio)
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
