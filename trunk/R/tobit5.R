tobit5c <- function( formula1, formula2, formula3,
                    data=sys.frame( sys.parent()), print.level=0, ...)
{
### Arvutab tobit5 t¸¸pi mudeli FIML abil (Amemiya 1985). See on
### tobit-5 mudel kujul:
### y1* = Z'g + u1
### y2* = X2'b + u2
### y3* = X3'b + u3
### S.o. mudel, kus koehvitsent on eri variantidel sama, normaalselt
### on ka X3 ja X2 samad.
###
### Parameetrid:
### formula1 - selektsioonimudel tavalisel R-i kujul
### formula2 - valik 2 mudel
### formula3 - valik 3 mudel
### NB! valikutel 2 ja 3 peab olema ¸hesugune mudeli dimensioon, nii
### et b saaks olla sama.
### Kıigi vektorite dimensioon peab olema sama, valiku 2 ja 3 vektorid peavad
### t‰ispikkusega. T‰ispikkusega vektori mittevaadeldavad juhud vıivad olla mis
### iganes. 
###  ...    - lisaparameetrid minimeerimisele
###  
### TULEMUSED:
### Tulemuseks on list, mille komponendid on:
### probit  - y1 probit koefitsentide maatriks
### H2S2    - Hekmani kahesammulise meetodiga leitud lahend teisele valikule.
###           See on list, mille komponendid on:
###  tulemused - OLS koefitsentide maatriks,
###  sigma     - j‰‰kliikmete dispersioon
### H2S3    - analoogiline lahend kolmandale valikule
### tulemused - ML keofitsentide ja veahinnangute maatriks
### minimeerimine - minimeerimisalgoritmi tulemused

  loglik <- function( beeta) {
                                        # Parameetrite lahti
                                        # kakkumine, p-d eelnevalt
                                        # defineeritud
    g <- beeta[pg]
    b <- beeta[pb]
    sigma2 <- beeta[psigma2]
    if( sigma2 <= 0) return( JAMA)
    roo2 <- beeta[proo2]
    if( ( roo2 < -1) || ( roo2 > 1)) return( JAMA)
    sigma3 <- beeta[psigma3]
    if( sigma3 <= 0) return( JAMA)
    roo3 <- beeta[proo3]
    if( ( roo3 < -1) || ( roo3 > 1)) return( JAMA)
                                        # Arvutamise kiirendamiseks
                                        # mıned suurused valmis:
    Z2g <- Z2%*%g
    Z3g <- Z3%*%g
    X2b2 <- X2%*%b
    X3b3 <- X3%*%b

    u2 <- y2 - X2b2
    u3 <- y3 - X3b3
    sqrt1r22 <- sqrt( 1 - roo2^2)
    sqrt1r32 <- sqrt( 1 - roo3^2)
    B2 <- -(Z2g + roo2/sigma2*u2)/sqrt1r22
    B3 <- (Z3g + roo3/sigma3*u3)/sqrt1r32
    fB2 <- dnorm( B2)
    FB2 <- pnorm( B2)
    fB3 <- dnorm( B3)
    FB3 <- pnorm( B3)
    lambda2 <- fB2/FB2
    lambda3 <- fB3/FB3
    CB2 <- as.vector( -(FB2*fB2*B2 + fB2*fB2)/FB2/FB2)
    CB3 <- as.vector( -(FB3*fB3*B3 + fB3*fB3)/FB3/FB3)
### log-laiklihuud    
    l2 <- sum(
              -log( sigma2) - 0.5*( u2/sigma2)^2 +
              pnorm( B2, log=TRUE))
    l3 <- sum(
              -log( sigma3) - 0.5*( u3/sigma3)^2 +
              pnorm( B3, log=TRUE))
    loglik <- l2 + l3 - T/2*log( 2*pi)
### N¸¸d gradient. S‰‰l on kokku 7 komponenti, neist 3
### vektorit. J‰rjekord j‰‰b samaks mis enne: g b_2, s_2, r_2, b_3,
### s_3, r_3
    l.g <- -t( Z2)%*%lambda2/sqrt1r22 + t( Z3)%*%lambda3/sqrt1r32
    l.b <- t( X2)%*%( lambda2*roo2/sigma2/sqrt1r22 + u2/sigma2^2) +
                                        # esimene komponent valik kahest
      t( X3)%*%( -lambda3*roo3/sigma3/sqrt1r32 + u3/sigma3^2)
                                        # teine komponent valik kolmest
    l.s2 <- sum( -1/sigma2 + u2^2/sigma2^3
                +lambda2*roo2/sigma2^2*u2/sqrt1r22)
    l.s3 <- sum( -1/sigma3 + u3^2/sigma3^3
                -lambda3*roo3/sigma3^2*u3/sqrt1r32)
    l.r2 <- -sum( lambda2*( u2/sigma2 + roo2*Z2g)/sqrt1r22^3)
    l.r3 <- sum( lambda3*( u3/sigma3 + roo3*Z3g)/sqrt1r32^3)
    gradient <- rbind( l.g, l.b, l.s2, l.r2, l.s3, l.r3)
### N¸¸d 28 Hessi maatriksi komponenti
    l.gg <- t( Z2) %*% ( Z2 * CB2)/sqrt1r22^2 +
      t( Z3) %*% ( Z3 * CB3)/sqrt1r32^2
    l.gb <- -t( Z2) %*%
      ( X2 * CB2)*roo2/sqrt1r22^2/sigma2 +
                                        # valik 2
        -t( Z3) %*%
          ( X3 * CB3)*roo3/sqrt1r32^2/sigma3
                                        # valik 3
    l.gs2 <- -roo2/sigma2^2/sqrt1r22^2*
      t( Z2) %*% ( CB2*u2)
    l.gr2 <- t( Z2) %*%
      ( CB2*( u2/sigma2 + roo2*Z2g)/sqrt1r22^4 -
       lambda2*roo2/sqrt1r22^3)
    l.gs3 <- -roo3/sigma3^2/sqrt1r32^2*
      t( Z3) %*% ( CB3*u3)
    l.gr3 <- t( Z3) %*%
      ( CB3*( u3/sigma3 + roo3*Z3g)/sqrt1r32^4 +
       lambda3*roo3/sqrt1r32^3)
    l.bb <- t( X2) %*%
      (X2 * ( (roo2/sqrt1r22)^2 * CB2 - 1))/sigma2^2 +
                                        # valik 2
        t( X3) %*%
          (X3 * ( (roo3/sqrt1r32)^2 * CB3 - 1))/sigma3^2
    l.bs2 <- t( X2) %*%
      ( CB2*roo2^2/sigma2^3*u2/sqrt1r22^2 -
       roo2/sigma2^2*lambda2/sqrt1r22 -
       2*u2/sigma2^3)
    l.br2 <- t( X2) %*%
      ( -CB2*( u2/sigma2 + roo2*Z2g)/sqrt1r22^4*roo2 +
       lambda2/sqrt1r22^3)/sigma2
    l.bs3 <- t( X3) %*%
      ( CB3*roo3^2/sigma3^3*u3/sqrt1r32^2 +
       roo3/sigma3^2*lambda3/sqrt1r32 - 2*u3/sigma3^3)
    l.br3 <- t( X3) %*%
      ( -CB3*( u3/sigma3 + roo3*Z3g)/sqrt1r32^4*roo3 -
       lambda3/sqrt1r32^3)/sigma3
    l.s2s2 <- sum(
                  1/sigma2^2
                  -3*u2*u2/sigma2^4
                  + u2*u2/sigma2^4 *roo2^2/sqrt1r22^2 *CB2
                  -2*lambda2* u2/sqrt1r22 *roo2/sigma2^3)
    l.s2r2 <- sum(
                  ( -CB2*roo2*(u2/sigma2 + roo2*Z2g)/sqrt1r22 +
                  lambda2)
                  *u2/sigma2^2)/sqrt1r22^3
    ## l.s2x3 on kıik nullid
    l.r2r2 <- sum(
                  CB2*( ( u2/sigma2 + roo2*Z2g)/sqrt1r22^3)^2
                  -lambda2*( Z2g*( 1 + 2*roo2^2) + 3*roo2*u2/sigma2) /
                    sqrt1r22^5
                  )
    ## l.r2x3 on kıik nullid
    l.s3s3 <- sum(
                  1/sigma3^2
                  -3*u3*u3/sigma3^4
                  +2*lambda3* u3/sqrt1r32 *roo3/sigma3^3
                  +roo3^2/sigma3^4 *u3*u3/sqrt1r32^2 *CB3)
    l.s3r3 <- -sum(
                  ( CB3*roo3*(u3/sigma3 + roo3*Z3g)/sqrt1r32 +
                  lambda3)
                  *u3/sigma3^2)/sqrt1r32^3
    l.r3r3 <- sum(
                  CB3*( ( u3/sigma3 + roo3*Z3g)/sqrt1r32^3)^2
                  + lambda3*( Z3g*( 1 + 2*roo3^2) + 3*roo3*u3/sigma3) /
                    sqrt1r32^5
                  )
    hess <- array( NA, c( df, df))
    hess[pg,pg] <- l.gg
    hess[pg,pb] <- l.gb; hess[pb,pg] <- t( l.gb)
    hess[pg,psigma2] <- l.gs2; hess[psigma2,pg] <- t( l.gs2)
    hess[pg,proo2] <- l.gr2; hess[proo2,pg] <- t( l.gr2)
    hess[pg,psigma3] <- l.gs3; hess[psigma3,pg] <- t( l.gs3)
    hess[pg,proo3] <- l.gr3; hess[proo3,pg] <- t( l.gr3)
    hess[pb,pb] <- l.bb
    hess[pb,psigma2] <- l.bs2; hess[psigma2,pb] <- t( l.bs2)
    hess[pb,proo2] <- l.br2; hess[proo2,pb] <- t( l.br2)
    hess[pb,psigma3] <- l.bs3; hess[psigma3,pb] <- t( l.bs3)
    hess[pb,proo3] <- l.br3; hess[proo3,pb] <- t( l.br3)
    hess[psigma2,psigma2] <- l.s2s2
    hess[psigma2,proo2] <- l.s2r2; hess[proo2,psigma2] <- l.s2r2
    hess[psigma2,psigma3] <- 0; hess[psigma3,psigma2] <- 0
    hess[psigma2,proo3] <- 0; hess[proo3,psigma2] <- 0
    hess[proo2,proo2] <- l.r2r2
    hess[proo2,psigma3] <- 0; hess[psigma3,proo2] <- 0
    hess[proo2,proo3] <- 0; hess[proo3,proo2] <- 0
    hess[psigma3,psigma3] <- l.s3s3
    hess[psigma3,proo3] <- l.s3r3; hess[proo3,psigma3] <- t( l.s3r3)
    hess[proo3,proo3] <- l.r3r3
    ll <- -loglik
    attr( ll, "gradient") <- -gradient
    attr( ll, "hessian") <- -hess
    ll
  }## loglik lıpp

  JAMA <- NA
  Z <- model.matrix( formula1, data=data)
  X2 <- model.matrix( formula2, data=data)
  X3 <- model.matrix( formula3, data=data)
                                        # X3 sisaldab
                                        # lisamuutujaid. lisaks paneme
                                        # talle hiljem lisakonstandi
  dfz <- ncol( Z)
  dfx <- ncol( X2)
  df <- dfz + dfx + 4
                                        # Parameetrite koguarv
  y1 <- model.response( model.frame( formula1, data=data))
  y2 <- model.response( model.frame( formula2, data=data))
  y3 <- model.response( model.frame( formula3, data=data))
  T <- length( y1)
  T2 <- length( y1[y1==0])
  T3 <- length( y1[y1==1])
                                        # N¸¸d leiame parameetrite
                                        # indeksid
  pg <- 1:dfz
  pb <- (dfz+1):(dfz+dfx)
                                        # yksiku mudeli koefitsendid
  psigma2 <- dfz+dfx+1
  proo2 <- dfz+dfx+2
  psigma3 <- dfz+dfx+3
  proo3 <- dfz+dfx+4
                                        # N¸¸d jagame andmed kahte
                                        # lehte, tolle j‰rgi, et kes
                                        # mis vıimaluse valib
  y2 <- y2[y1==0]
  y3 <- y3[y1==1]
  X2 <- X2[y1==0,]
  X3 <- X3[y1==1,]
                                        # NB! X2 ja X3 siin eri pikkusega
  Z2 <- Z[y1==0,]
  Z3 <- Z[y1==1,]
  cat( "Valik 2:", T2, "korda; valik 3:", T3, "korda\n")
  b0 <- numeric()
### Algv‰‰rtused. Probit osa algv‰‰rtused probitiga
  probit <- probit( formula1, data=data)
  if( print.level > 0) {
    cat( "Probit osa algv‰‰rtused (probit mudel):\n")
    print( probit$tulemused)
  }
  gamma.0 <- probit$tulemused[,1]
  b0[pg] <- gamma.0
                                        # kirjutame gamma algv‰‰rtuse
                                        # parameetrite vektorisse
  ## mılemad muutujad Heckmani kahesammulise meetodiga
  lambda.2 <- dnorm( -Z2%*%gamma.0)/pnorm( -Z2%*%gamma.0)
  lambda.3 <- dnorm( Z3%*%gamma.0)/pnorm( Z3%*%gamma.0)
  delta.2 <- mean( lambda.2^2 - Z2%*%gamma.0 *lambda.2)
                                        # delta on OLS dispersiooni kordaja
  delta.3 <- mean( lambda.3^2 + Z3%*%gamma.0 *lambda.3)
  kahesammuline.2 <- summary( lm( y2 ~ -1+X2+lambda.2))
                                        # X2 sisaldab juba konstanti
  kahesammuline.3 <- summary( lm( y3 ~ -1+X3+lambda.3))
  b0[pb] <- kahesammuline.2$coefficients[1:dfx,1]
                                        # yksiku mudeli koefitsendid
  se.2 <- kahesammuline.2$sigma
  se.3 <- kahesammuline.3$sigma
                                        # OLS dispersioonid
  bl.2 <- kahesammuline.2$coefficients[dfx+1,1]
  bl.3 <- kahesammuline.3$coefficients[dfx+1,1]
                                        # OLS lambda kordajad
  sigma.20 <- sqrt( se.2^2 + ( bl.2*delta.2)^2)
  sigma.30 <- sqrt( se.3^2 + ( bl.3*delta.3)^2)
  b0[psigma2] <- sigma.20
  b0[psigma3] <- sigma.30
  roo2.0 <- -bl.2/sigma.20
  roo3.0 <- bl.3/sigma.30
  if( roo2.0 <= -1) roo2.0 <- -0.99
  if( roo3.0 <= -1) roo3.0 <- -0.99
  if( roo2.0 >= 1) roo2.0 <- 0.99
  if( roo3.0 >= 1) roo3.0 <- 0.99
  b0[proo2] <- roo2.0
  b0[proo3] <- roo3.0
  if( print.level > 0) {
    cat( "beeta algv‰‰rtused (Heckmani kahesammuline mudel):\n")
    cat( b0[pb], "\n")
    cat( "sigma2, roo2; sigma3, roo3:\n")
    cat( b0[psigma2], " ", b0[proo2], "; ", b0[psigma3], " ", b0[proo3], "\n")
  }
  a <- minNR( loglik, b0, print.level=print.level, ...)
  if( a$code > 3) {
    cat( "Jama: es leia lahendit -", a$kommentaar, "\n")
    stdh <- 0
  }
  else {
    stdh <- sqrt( diag( solve( a$hessian)))
                                        # miinust pole, sest loglik
                                        # antakse nagunii miinusm‰rgiga
  }
  t <- a$estimate/stdh
  p <- 2*pnorm( -abs( t))
  tulemused <<- cbind( a$estimate, stdh, t, p)

  Xmuutujate.nimed <- dimnames( X3)[[2]]
  Zmuutujate.nimed <- dimnames( Z)[[2]]
  Xmuutujate.nimed[1] <- "konst"
  Zmuutujate.nimed[1] <- "konst"
  dimnames( tulemused) <- list( as.list( c( Zmuutujate.nimed, 
                                           Xmuutujate.nimed,
                                           "s2", "r2", "s3", "r3")),
                               list( "koefitsent", "stdh", "t", "P(|b|>t)"))
  tobit5.tulemus <- list( probit=probit$tulemused,
                         H2S2=list(
                            tulemused=kahesammuline.2$coefficients,
                            sigma=kahesammuline.2$sigma),
                          H2S3=list(
                            tulemused=kahesammuline.3$coefficients,
                            sigma=kahesammuline.3$sigma),
                         tulemused=tulemused, minimeerimine=a)
}
