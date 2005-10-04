summary.maximisation <- function(object, hessian=FALSE, unsucc.step=FALSE,
   ... ) {
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
   estimate <- cbind(object$estimate, object$gradient)
   if(hessian) {
      H <- object$hessian
   }
   else {
      H <- NULL
   }
   summary <- list(type=object$type,
                  iterations=object$iterations,
                  code=object$code,
                  message=object$message,
                   unsucc.step=unsucc.step,
                  estimate=estimate,
                   hessian=H)
   class(summary) <- "summary.maximisation"
   summary
}
