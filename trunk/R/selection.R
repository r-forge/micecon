selection <- function(selection, outcome,
                      data=sys.frame(sys.parent()),
                      method="ml",
                      init=NULL, print.level=0,
                      ySelection=FALSE, xSelection=FALSE,
                      yOutcome=FALSE, xOutcome=FALSE,
                      model=FALSE,
                      ...) {
   ## Heckman-style sample-selection models
      ## --- the main program ---
   ## First the consistency checks
   if( class( selection ) != "formula" ) {
      stop( "argument 'outcome' must be a formula" )
   }
   if( length( formula ) != 3 ) {
      stop( "argument 'outcome' must be a 2-sided formula" )
   }
   if( "probit" %in% substr( all.vars( outcome ), 1, 6 ) ) {
      stop( "argument 'outcome' may not include variable names",
                  " starting with 'probit'" )
   }
   if( class( selection ) != "formula" ) {
      stop( "argument 'selection' must be a formula" )
   }
   if( length( selection ) != 3 ) {
      stop( "argument 'selection' must be a 2-sided formula" )
   }
   if( "probit" %in% substr( all.vars( selection ), 1, 6 ) ) {
      stop( "argument 'selection' may not include a variable",
                  " names starting with 'probit'" )
   }
   probitEndogenous <- model.frame( selection, data = data )[ , 1 ]
   probitLevels <- levels( as.factor( probitEndogenous ) )
   if( length( probitLevels ) != 2 ) {
      stop( "the left hand side of 'selection' has to contain",
         " exactly two levels (e.g. FALSE and TRUE)" )
   }
   data$probitDummy <- probitEndogenous == probitLevels[ 2 ]
   ## now check whether two-step method is needed: either for final estimate or
   ## initial parameters
   if(method == "2step" | is.null(init)) {
      twoStep <- heckit(selection, outcome, data)
      if(method == "2step") {
         return(twoStep)
      }
   }
   ## Now extract model frames etc
   ## YS (selection equation)
   cl <- match.call()
   mf <- match.call(expand.dots = FALSE)
   m <- match(c("selection", "data", "subset", "weights", "na.action",
                "offset"), names(mf), 0)
   mf1 <- mf[c(1, m)]
   mf1$drop.unused.levels <- TRUE
   mf1[[1]] <- as.name("model.frame")
   names(mf1)[2] <- "formula"
                                        # model.frame requires the parameter to
                                        # be 'formula'
   mf1 <- eval(mf1, parent.frame())
   mt1 <- attr(mf1, "terms")
   XS <- model.matrix(mt1, mf1)
   YS <- model.response(mf1, "numeric")
   ## YO (outcome equation)
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
   ## now fit the model
   estimation <- tobit2.fit(YS, XS, YO, XO, init)
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
               terms1=mt1,
               terms2=mt2,
               y1=switch(y1, "1"=list(Y1), "0"=NULL),
               z=switch(z, "1"=list(X), "0"=NULL),
               y2=switch(y2, "1"=list(Y2), "0"=NULL),
               x=switch(x, "1"=list(X), "0"=NULL),
               model=switch(model, "1"=list(selection=mf1, formula=mf2), "0"=NULL))
   class(result) <- c("sampleSelection", "tobit2", class(estimation))
   return(result)
}

   
