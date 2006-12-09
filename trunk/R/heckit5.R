heckit5 <- function(selection, outcome1, outcome2,
                    data=sys.frame(sys.parent()),
                    ySelection=FALSE, xSelection=FALSE,
                    yOutcome=FALSE, xOutcome=FALSE,
                    model=FALSE,
                    print.level=0,
                    ...) {
   ## Do a few sanity checks...
   if( class( selection ) != "formula" ) {
      stop( "argument 'selection' must be a formula" )
   }
   if( length( selection ) != 3 ) {
      stop( "argument 'selection' must be a 2-sided formula" )
   }
   if("formula" %in% class( outcome1 )) {
      if( length( outcome1 ) != 3 ) {
         stop( "argument 'outcome1' must be a 2-sided formula" )
      }
   }
   else
       stop("argument 'outcome1' must be a formula")
   if("formula" %in% class( outcome2 )) {
      if( length( outcome2 ) != 3 ) {
         stop( "argument 'outcome2' must be a 2-sided formula" )
      }
   }
   else
       stop("argument 'outcome2' must be a formula")
   probitEndogenous <- model.frame( selection, data = data )[ , 1 ]
   probitLevels <- levels( as.factor( probitEndogenous ) )
   if( length( probitLevels ) != 2 ) {
      stop( "the left hand side of 'selection' has to contain",
         " exactly two levels (e.g. FALSE and TRUE)" )
   }
   ## Now extract model frames etc.  We only need model information for the selection equation
   cl <- match.call()
   mf <- match.call(expand.dots = FALSE)
   m <- match(c("selection", "data", "subset", "weights", "na.action",
                "offset"), names(mf), 0)
   mfS <- mf[c(1, m)]
   mfS$drop.unused.levels <- TRUE
   mfS[[1]] <- as.name("model.frame")
   names(mfS)[2] <- "formula"
                                        # model.frame requires the parameter to
                                        # be 'formula'
   mfS <- eval(mfS, parent.frame())
   mtS <- attr(mfS, "terms")
   XS <- model.matrix(mtS, mfS)
   YS <- factor(model.response(mfS, "numeric"))
   XS1 <- XS[YS == levels(YS)[1],]
   XS2 <- XS[YS == levels(YS)[2],]
   ## and run the model
   probit <- probit(selection)
   if( print.level > 0) {
      cat("The probit part of the model:\n")
      print(summary(probit))
   }
   gamma <- coefficients(probit)
   ## regression parts by Heckman's two-step method
   lambda1 <- dnorm( -XS1%*%gamma)/pnorm( -XS1%*%gamma)
   lambda2 <- dnorm( XS2%*%gamma)/pnorm( XS2%*%gamma)
                                        # lambdas are inverse Mills ratios
   twoStep1 <- lm(update(outcome1, ~ . + lambda1), data=data)
   twoStep2 <- lm(update(outcome2, ~ . + lambda2), data=data)
   se1 <- twoStep1$sigma
   se2 <- twoStep2$sigma
                                        # residual variance
   delta1 <- mean( lambda1^2 - XS1%*%gamma *lambda1)
   delta2 <- mean( lambda2^2 + XS2%*%gamma *lambda2)
   betaL1 <- coef(twoStep1)["lambda1"]
   betaL2 <- coef(twoStep2)["lambda2"]
   sigma1 <- sqrt( se1^2 + ( betaL1*delta1)^2)
   sigma2 <- sqrt( se2^2 + ( betaL2*delta2)^2)
   rho1 <- -betaL1/sigma1
   rho2 <- betaL2/sigma2
   if( rho1 <= -1) rho1 <- -0.99
   if( rho2 <= -1) rho2 <- -0.99
   if( rho1 >= 1) rho1 <- 0.99
   if( rho2 >= 1) rho2 <- 0.99
   result <- list(probit=probit,
                  twoStep1=twoStep1,
                  rho1=rho1,
                  sigma1=sigma1,
                  twoStep2=twoStep2,
                  rho2=rho2,
                  sigma2=sigma2)
}
