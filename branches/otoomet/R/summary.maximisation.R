summary.maximisation <- function( object, hessian=FALSE, unsucc.step=FALSE, ... ) {
   ### The object of class "maximisation" should include following components:
   ### maximum    : function value at optimum
   ### estimate   : estimated parameter values at optimum
   ### gradient   :           gradient at optimum
   ### hessian    :           hessian
   ### code       : code of convergence
   ### message    : message, description of the code
   ### last.step  : information about last step, if unsuccessful
   ### iterations : number of iterations
   ### type       : type of optimisation
   ###
   NParam <- length(object$estimate)
   if(!is.null(object$acivePar)) {
      activePar <- object$activePar
   } else {
      activePar <- rep(TRUE, NParam)
   }
   if(object$code == 3 & unsucc.step) {
      a <- cbind(object$last.step$theta0, object$last.step$theta1)
      dimnames(a) <- list(parameter=object$names,
                        c("current par", "new par"))
      unsucc.step <- list(value=object$last.step$f0,
                        parameters=a)
   } else {
      unsucc.step <- NULL
   }
   if(object$code < 100) {
      if(min(abs(eigen(object$hessian[activePar,activePar],
            symmetric=TRUE, only.values=TRUE)$values)) > 1e-6) {
         varcovar <- matrix(0, NParam, NParam)
         varcovar[activePar,activePar] <-
            solve(-object$hessian[activePar,activePar])
         hdiag <- diag(varcovar)
         if(any(hdiag < 0)) {
            warning("Diagonal of variance-covariance matrix not positive!\n")
         }
         stdh <- sqrt(hdiag)
         t <- object$estimate/stdh
         p <- 2*pnorm( -abs( t))
      } else {
         stdh <- 0
         t <- 0
         p <- 0
      }
      if(!is.null(object$gradient)) {
         results <- cbind(object$estimate, stdh, t, p,
                       as.numeric(object$gradient))
         dimnames(results) <- list(parameter=names(object$estimate),
                                list("coeff.", "stdd", "t",
                                     "P(|b|>t)", "gradient"))
      } else {
         results <- cbind(object$estimate, stdh, t, p, 1 - as.integer(activeParam))
         dimnames(results) <- list(parameter=object$names,
                                list("coeff.", "stdd", "t",
                                     "P(|b|>t)", "constant"))
      }
      estimate <- list(value=object$maximum,
                     results=results)
      if(hessian) {
         estimate <- c(estimate, list(hessian=object$hessian))
      }
   } else {
      estimate <- NULL
   }
   summary <- list(type=object$type,
                  iterations=object$iterations,
                  code=object$code,
                  message=object$message,
                  unsucc.step=unsucc.step,
                  estimate=estimate)
   class(summary) <- "summary.maximisation"
   summary
}
