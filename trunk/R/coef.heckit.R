coef.heckit <- function(object, part="outcome", ...) {
   if(part=="outcome") {
      b <- object$coefficients[c(object$param$index$betaO, object$param$index$invMillsRatio)]
   }
   else if(part=="full") {
      b <- object$coefficients
   }
   else
       stop("'part' must be either 'outcome' or 'full'")
   b
}
