heckit5 <- function(selection, outcome1, outcome2,
                    data=sys.frame(sys.parent()),
                    ys=FALSE, yo=FALSE,
                    xs=FALSE, xo=FALSE,
                    mfs=FALSE, mfo=FALSE,
                    print.level=0) {
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
   i1 <- YS == levels(YS)[1]
   i2 <- YS == levels(YS)[2]
   XS1 <- XS[i1,]
   XS2 <- XS[i2,]
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
   colnames(lambda1) <- colnames(lambda2) <- "invMillsRatio"
                                        # lambdas are inverse Mills ratios
   ## We need here to run OLS on outcome1 expression and lambda1.  Outcome1 must be evaluated in the
   ## parent frame, lambda1 here -> we must evaluate it here
   m <- match(c("outcome1", "data", "subset", "weights", "na.action",
                "offset"), names(mf), 0)
   mf1 <- mf[c(1, m)]
   mf1$drop.unused.levels <- TRUE
   mf1[[1]] <- as.name("model.frame")
   names(mf1)[2] <- "formula"
   mf1 <- eval(mf1, parent.frame())
   mt1 <- attr(mf1, "terms")
   X1 <- model.matrix(mt1, mf1)[i1,]
   Y1 <- model.response(mf1, "numeric")[i1]
   X1 <- cbind(X1, lambda1)
                                        # lambda1 is a matrix -- we need to remove the dim in order to 
   m <- match(c("outcome2", "data", "subset", "weights", "na.action",
                "offset"), names(mf), 0)
   mf2 <- mf[c(1, m)]
   mf2$drop.unused.levels <- TRUE
   mf2[[1]] <- as.name("model.frame")
   names(mf2)[2] <- "formula"
   mf2 <- eval(mf2, parent.frame())
   mt2 <- attr(mf2, "terms")
   X2 <- model.matrix(mt2, mf2)[i2,]
   Y2 <- model.response(mf2, "numeric")[i2]
   X2 <- cbind(X2, lambda2)
                                        #
   twoStep1 <- lm(Y1 ~ -1 + X1)
   twoStep2 <- lm(Y2 ~ -1 + X2)
                                        # X includes the constant
   se1 <- summary(twoStep1)$sigma
   se2 <- summary(twoStep2)$sigma
                                        # residual variance
   delta1 <- mean( lambda1^2 - XS1%*%gamma *lambda1)
   delta2 <- mean( lambda2^2 + XS2%*%gamma *lambda2)
   betaL1 <- coef(twoStep1)["X1invMillsRatio"]
   betaL2 <- coef(twoStep2)["X2invMillsRatio"]
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
                  sigma2=sigma2,
                  termsS=mtS,
                  termsO=list(mt1, mt2),
                  ys=switch(as.character(ys), "TRUE"=YS, "FALSE"=NULL),
                  xs=switch(as.character(xs), "TRUE"=XS, "FALSE"=NULL),
                  yo=switch(as.character(yo), "TRUE"=list(YO1, YO2), "FALSE"=NULL),
                  xo=switch(as.character(xo), "TRUE"=list(XO1, XO2), "FALSE"=NULL),
                  mfs=switch(as.character(mfs), "TRUE"=list(mfS), "FALSE"=NULL),
                  mfo=switch(as.character(mfs), "TRUE"=list(mf1, mf2), "FALSE"=NULL)
                  )
   class(result) <- c("heckit5", "heckit", class(result))
   invisible(result)
}
