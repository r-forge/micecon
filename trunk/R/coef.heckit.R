coef.heckit <- function(object, part="outcome", ...) {
   if(part=="outcome") {
      b <- object$coefficients[c(object$param$index$betaO, object$param$index$Mills)]
   }
   else if(part=="full") {
      b <- object$coefficients
   }
   else
       stop("'part' must be either 'outcome' or 'full'")
   b
}

coef.heckit5 <- function(object, part="outcome", ...) {
   if(part=="outcome") {
      b <- object$coefficients[c(object$param$index$betaO1, object$param$index$invMillsRatio1,
                                 object$param$index$betaO2, object$param$index$invMillsRatio2)]
   }
   else if(part=="full") {
      b <- object$coefficients
   }
   else
       stop("'part' must be either 'outcome' or 'full'")
   b
}
