selection <- function(selection, outcome,
                      data=sys.frame(sys.parent()),
                      method="ml",
                      init=NULL,
                      ySelection=FALSE, xSelection=FALSE,
                      yOutcome=FALSE, xOutcome=FALSE,
                      model=FALSE,
                      print.level=0,
                      ...) {
   ## Heckman-style sample-selection models
   ## selection:   formula
   ##              LHS: must be convertable to two-level factor (e.g. 0-1, 1-2, "A"-"B")
   ##              RHS: ordinary formula as in lm()
   ## outcome:     formula, or a list of two formulas.
   ##              If the outcome contains one formula, it is
   ##              assumed that we observe the formula only if selection equals to the first level
   ##              of selection (e.g. "0", "1" or "A" in the examples above).
   ##              If the outcome contains two formulas, we assume that the first one is observed if
   ##              selection equals to the first level, otherwise the second formula.
   ## First the consistency checks
   type <- 0
   if( class( selection ) != "formula" ) {
      stop( "argument 'selection' must be a formula" )
   }
   if( length( selection ) != 3 ) {
      stop( "argument 'selection' must be a 2-sided formula" )
   }
                                        #
   if("formula" %in% class( outcome )) {
      if( length( outcome ) != 3 ) {
         stop( "argument 'outcome' must be a 2-sided formula" )
      }
      type <- 2
   }
   else if("list" %in% class(outcome)) {
      if(length(outcome) == 1) {
         outcome <- outcome[[1]]
         type <- 2
      }
      else if(length(outcome) == 2) {
         if("formula" %in% class(outcome[[1]])) {
            if( length( outcome[[1]] ) != 3 ) {
               stop( "argument 'outcome[[1]]' must be a 2-sided formula" )
            }
         }
         else
             stop( "argument 'outcome[[1]]' must be either a formula or a list of two formulas" )
         if("formula" %in% class(outcome[[2]])) {
            if( length( outcome[[2]] ) != 3 ) {
               stop( "argument 'outcome[[2]]' must be a 2-sided formula" )
            }
         }
         else
             stop( "argument 'outcome[[2]]' must be either a formula or a list of two formulas" )
         type <- 5
      }
      else
          stop("argument 'outcome' must contain 1 or 2 components")
   }
   else
       stop( "argument 'selection' must be either a formula or a list of two formulas" )
   if(print.level > 0)
       cat("Tobit", type, "model\n")
   probitEndogenous <- model.frame( selection, data = data )[ , 1 ]
   probitLevels <- levels( as.factor( probitEndogenous ) )
   if( length( probitLevels ) != 2 ) {
      stop( "the left hand side of 'selection' has to contain",
         " exactly two levels (e.g. FALSE and TRUE)" )
   }
   data$probitDummy <- probitEndogenous == probitLevels[ 2 ]
   ## now check whether two-step method was requested
   if(method == "2step") {
      if(type == 2)
          twoStep <- heckit2(selection, outcome, data)
      else if(type == 5)
          twoStep <- heckit5(selection, outcome, data)
      else
          stop("unknown type")
      return(twoStep)
   }
   ## Now extract model frames etc
   ## YS (selection equation)
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
   YS <- model.response(mfS, "numeric")
   ## YO (outcome equation)
   if(type == 2) {
      m <- match(c("outcome", "data", "subset", "weights", "na.action",
                   "offset"), names(mf), 0)
      mf2 <- mf[c(1, m)]
      mf2$drop.unused.levels <- TRUE
      mf2[[1]] <- as.name("model.frame")
      names(mf2)[2] <- "formula"
      mf2 <- eval(mf2, parent.frame())
                                        # Note: if unobserved variables are
                                        # marked as NA, eval returns a
                                        # subframe of visible variables only.
                                        # We have to check it later
      mt2 <- attr(mf2, "terms")
      XO <- model.matrix(mt2, mf2)
      YO <- model.response(mf2, "numeric")
      estimation <- tobit2fit(YS, XS, YO, XO, init)
   }
   if(type == 5) {
      oArg <- match("outcome", names(mf), 0)
                                        # find the outcome argument
      ocome <- as.list(mf[[oArg]])
      formula1 <- ocome[[2]]
      formula2 <- ocome[[3]]
                                        # Now we have extracted both formulas
      m <- match(c("outcome", "data", "subset", "weights", "na.action",
                   "offset"), names(mf), 0)
      ## replace the outcome list by the first equation and evaluate it
      mf[[oArg]] <- formula1
      mf1 <- mf[c(1, m)]
      mf1$drop.unused.levels <- TRUE
      mf1[[1]] <- as.name("model.frame")
                                        # eval it as model frame
      names(mf1)[2] <- "formula"
      mf1 <- eval(mf1, parent.frame())
                                        # Note: if unobserved variables are
                                        # marked as NA, eval returns a
                                        # subframe of visible variables only.
                                        # We have to check it later
      mt1 <- attr(mf1, "terms")
      XO1 <- model.matrix(mt1, mf1)
      YO1 <- model.response(mf1, "numeric")
      ## repeat all the stuff with second equation
      mf[[oArg]] <- formula2
      mf2 <- mf[c(1, m)]
      mf2$drop.unused.levels <- TRUE
      mf2[[1]] <- as.name("model.frame")
                                        # eval it as model frame
      names(mf2)[2] <- "formula"
      mf2 <- eval(mf2, parent.frame())
                                        # Note: if unobserved variables are
                                        # marked as NA, eval returns a
                                        # subframe of visible variables only.
                                        # We have to check it later
      mt2 <- attr(mf2, "terms")
      XO2 <- model.matrix(mt2, mf2)
      YO2 <- model.response(mf2, "numeric")
      ## indices in for the parameter vector.  These are returned in order to provide the user a way
      ## to extract certain components from the coefficients
      NXS <- ncol(XS)
      NXO1 <- ncol(XO1)
      NXO2 <- ncol(XO2)
      igamma <- 1:NXS
      ibeta1 <- seq(tail(igamma, 1)+1, length=NXO1)
      isigma1 <- tail(ibeta1, 1) + 1
      irho1 <- tail(isigma1, 1) + 1
      ibeta2 <- seq(tail(irho1, 1) + 1, length=NXO2)
      isigma2 <- tail(ibeta2, 1) + 1
      irho2 <- tail(isigma2, 1) + 1
      NParam <- irho2
      if(is.null(init)) {
         init <- numeric(NParam)
         heckit <- heckit5(selection, as.formula(formula1), as.formula(formula2),
                           data, print.level)
         init[igamma] <- coef(heckit$probit)
         b1 <- coef(heckit$twoStep1)
         init[ibeta1] <- b1[names(b1) != "invMillsRatio"]
         init[isigma1] <- heckit$sigma1
         init[irho1] <- heckit$rho1
         b2 <- coef(heckit$twoStep2)
         init[ibeta2] <- b2[names(b2) != "invMillsRatio"]
         init[isigma2] <- heckit$sigma2
         init[irho2] <- heckit$rho2
      }
      estimation <- tobit5fit(YS, XS, YO1, XO1, YO2, XO2, init)
   }
   ## now fit the model
   result <- c(estimation,
               twoStep=list(twoStep),
               NParam=NParam,
               NObs=NObs,
               N1=N1,
               N2=N2,
               NZ=NZ,
               NX=NX,
               df=NObs - NParam,
               call=cl,
               terms1=mtS,
               terms2=mt2,
               y1=switch(y1, "1"=list(Y1), "0"=NULL),
               z=switch(z, "1"=list(X), "0"=NULL),
               y2=switch(y2, "1"=list(Y2), "0"=NULL),
               x=switch(x, "1"=list(X), "0"=NULL),
               model=switch(model, "1"=list(selection=mfS, formula=mf2), "0"=NULL))
   class(result) <- c("sampleSelection", "tobit2", class(estimation))
   return(result)
}

   
