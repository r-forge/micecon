summary.maximisation <- function(result, hessian=FALSE, unsucc.step=FALSE) {
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
   NParam <- length(result$estimate)
   if(!is.null(result$acivePar)) {
      activePar <- result$activePar
   } else {
      activePar <- rep(TRUE, NParam)
   }
   if(result$code == 3 & unsucc.step) {
      a <- cbind(result$last.step$theta0, result$last.step$theta1)
      dimnames(a) <- list(parameter=result$names,
                        c("current par", "new par"))
      unsucc.step <- list(value=result$last.step$f0,
                        parameters=a)
   } else {
      unsucc.step <- NULL
   }
   if(result$code < 100) {
      if(min(abs(eigen(result$hessian[activePar,activePar],
            symmetric=TRUE, only.values=TRUE)$values)) > 1e-6) {
         varcovar <- matrix(0, NParam, NParam)
         varcovar[activePar,activePar] <-
            solve(-result$hessian[activePar,activePar])
         hdiag <- diag(varcovar)
         if(any(hdiag < 0)) {
            warning("Diagonal of variance-covariance matrix not positive!\n")
         }
         stdh <- sqrt(hdiag)
         t <- result$estimate/stdh
         p <- 2*pnorm( -abs( t))
      } else {
         stdh <- 0
         t <- 0
         p <- 0
      }
      if(!is.null(result$gradient)) {
         results <- cbind(result$estimate, stdh, t, p,
                       as.numeric(result$gradient))
         dimnames(results) <- list(parameter=names(result$estimate),
                                list("coeff.", "stdd", "t",
                                     "P(|b|>t)", "gradient"))
      } else {
         results <- cbind(result$estimate, stdh, t, p, 1 - as.integer(activeParam))
         dimnames(results) <- list(parameter=result$names,
                                list("coeff.", "stdd", "t",
                                     "P(|b|>t)", "constant"))
      }
      estimate <- list(value=result$maximum,
                     results=results)
      if(hessian) {
         estimate <- c(estimate, list(hessian=result$hessian))
      }
   } else {
      estimate <- NULL
   }
   summary <- list(type=result$type,
                  iterations=result$iterations,
                  code=result$code,
                  message=result$message,
                  unsucc.step=unsucc.step,
                  estimate=estimate)
   class(summary) <- "summary.maximisation"
   summary
}
