probit <- function( formula, b0=NULL, data=sys.frame(sys.parent()),
                   print.level=0, ...) {
### Probit binary choice model
### formula: model formula, response must be logical or 0/1 variable
### b0       initial value of the parameters.
### marginal - kas anda ka piirefektid
### ... - argumendid minimeerimisalgoritmile
###
### vastuseks annab listi ,probit.tulemus', mis sisaldab:
###  $tulemused: Nx4 maatriks, mis sisaldab koeftisente,
###   standardhälbeid, t-väärtuse ja tõenäosusi
###  $LRT - laiklihuud reitiõu test tolle kohta, et kogu probit mudel
###   on jama
###  $minimeerimine - minimeerimisalgoritmi tulemus
###  $marginal - piirefektid (juhul kui)
###
    loglik <- function( beta) {
        xb0 <- as.vector( x0 %*% beta)
        xb1 <- as.vector( x1 %*% beta)
        logF0 <- pnorm( xb0, log=TRUE)
        logF1 <- pnorm( xb1, log=TRUE)
        loglik <- sum( logF1) + sum( logF0)
    }
    gradlik <- function(beta) {
        xb0 <- as.vector( x0 %*% beta)
        xb1 <- as.vector( x1 %*% beta)
        F0 <- pnorm( xb0)
        F1 <- pnorm( xb1)
        f0 <- dnorm( xb0)
        f1 <- dnorm( xb1)
        yF0 <- pnorm( xb0, lower.tail=FALSE)
        loglik1 <- t( x1) %*% ( f1/F1) - t( x0) %*% ( f0/yF0)
    }
    hesslik <- function(beta) {
        xb0 <- as.vector( x0 %*% beta)
        xb1 <- as.vector( x1 %*% beta)
        F0 <- pnorm( xb0)
        F1 <- pnorm( xb1)
        f0 <- dnorm( xb0)
        f1 <- dnorm( xb1)
        yF0 <- pnorm( xb0, lower.tail=FALSE)
        loglik2 <- t( x1) %*% ( x1 * ( -f1*xb1*F1 - f1*f1)/F1/F1) -
            t( x0) %*% ( x0 * ( -f0*xb0*yF0 + f0*f0)/yF0/yF0)
                                        # note that df/db' = -f
                                        # (x'b) x'
    }
    x <- model.matrix( formula, data=data)
    y <- as.numeric(model.response( model.frame( formula, data=data)))
    df <- ncol( x)
                                        # df on X koos konstandiga
    T <- length( y)
    N1 <- length( y[y==1])
    x0 <- x[y==0,]
    x1 <- x[y==1,]
                                        # partial data, for speed
    if(is.null(b0))
        b0 <- rep( 0, df)
    ## Main estimation
    estimation <- maxNR(loglik, gradlik, hesslik, b0, print.level, ...)
    ## Now names for the variables
    muutujate.nimed <- dimnames(x)[[2]]
    muutujate.nimed[1] <- "konst"
    names(estimation$estimate) <- muutujate.nimed
    ## Likelihood ratio test: H0 -- all the coefficients, except intercept
    ## are zeros.  ML estimate for the intercept is qnorm(N1/T)
    ll.bar <- loglik(c(qnorm(N1/T), rep(0, df-1)))
    LRT <- 2*(estimation$maximum - ll.bar)
### Lõppvastus
    result <- list(results=estimation,
                   LRT=list(LRT=LRT, df=df-1))
                                        # there are df-1 constraints
    class(result) <- "probit"
    result
}

