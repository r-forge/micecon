piecewise.discrete <- function(t, censored, X=NULL, cm, b0=NULL,
                               normalise=FALSE,
                               ... ){
   ### t - measured time
   ### c - indicatiors for censoring 0 - uncensored, 1 - censored
   ### X - covariates w/o constant!
   ### cm - start times for intervals (the first is 0 not to include)
   ### b0 - initial parameters (optional)
   ### normalise - whether to require vl < 0 in solution
   ##
   ti <- function(t, cm) {
      ## time indeks - which interval the time t falls to (it falls into
      ## interval cm[ti] .. cm[ti+1].
      ## t - Nx1 vector
      ## cm - Mx1 vector
      ti <- rep(1, length(t))
      for( i in seq(cm)) {
         ti[t > cm[i]] <- i
      }
      invisible(ti)
   }
   ##
   loglik <- function(beta) {
      ### parameter vector beta:
      ### vl   - value of low hazard (normalised to exp(vl) and vh = 1)
      ### pl   - probability of low hazard (normalised to plogis(pl))
      ### lambda - baseline hazard values for time pieces (exp(lambda))
      ### beta   - coefficients for individual-specific hazard (exp(beta%*%X))
      vl <- exp(beta[ivl])
      # vl <- beta[ivl]
      # if(vl < 0) return(NA)
      pl <- plogis(beta[ipl])
      # pl <- beta[ipl]
      if(pl < 0 | pl > 1) {
         return(NA)
      }
      ph <- 1 - pl
      lambda <- beta[ilambda]
      gamma <- beta[igamma]
      if(!is.null(X)) {
         bX <- X %*% gamma
         mu <- exp(bX)
      } else {
         bX <- 0
         mu <- rep(1, N)
      }
      mu <- as.vector(exp(bX))
      surv <- matrix(c(diff(cm),0), N, M, byrow=TRUE)*cbind(Di, 0)
      surv[cbind(1:N,d)] <- t - cm[d]
                                          # columns: time survived in particular
                                          # pieces.
      lambda.pieces <- matrix(exp(lambda), N, M, byrow=TRUE)
                                          # columns: corresponding baseline
                                          # hazard
      hazard.pieces <- surv*lambda.pieces
      z <- as.vector(mu*(hazard.pieces %*% rep(1,M)))
                                          # survived time in the pieces *
                                          # corresponding baseline hazard, summed
                                          # over the pieces
      z.l <- as.vector(mu*vl*((hazard.pieces) %*% rep(1, M)))
                                          # the same for low value of v
      exit.hazard <- ifelse(NC, mu*exp(lambda[d]), 1)
                                          # hazard connected with exit from the
                                          # state.  Present only if uncensored
      exit.hazard.l <- ifelse(NC, vl*mu*exp(lambda[d]), 1)
      di.hazard <- di*matrix(exp(lambda), N, M, byrow=TRUE)
      ##
      likelihood.l <- as.vector(exit.hazard.l*exp(-z.l))
      likelihood.h <- as.vector(exit.hazard*exp(-z))
      likelihood <- pl*likelihood.l + (1-pl)*likelihood.h
      LR.l <- likelihood.l/likelihood
      LR.h <- likelihood.h/likelihood
                                          # NB! those were vectors
      loglik <- sum(log(likelihood))

      ## now gradient
      grad <- numeric(K)
      grad[ivl] <- sum(pl*LR.l*(NC/vl - z)) * vl
                                          # last term is associated with
                                          # normalization
      grad.ipl <- sum(LR.l - LR.h)
                                          # gradient in internal pl
      grad[ipl] <- grad.ipl* pl*(1 - pl)
                                          # and normalization for external one
      dLl.dlL <- LR.l*(NC*di - vl*mu*hazard.pieces)
                                          # 1/likelihood * dlikelihood.l/dlambda
      dLh.dlL <- LR.h*(NC*di - mu*hazard.pieces)
      vgrad.lambda <- pl*dLl.dlL + ph*dLh.dlL
      grad[ilambda] <- t(vgrad.lambda) %*% rep(1,N)

      if(!is.null(X)) {
         vgrad.gamma <- pl*LR.l*(NC - z.l) + ph*LR.h*(NC - z)
         grad[igamma] <- t(X) %*% vgrad.gamma
      }

      ## now hessian
      hess <- matrix(NA, K, K)
      hess[ivl,ivl] <- pl*sum(LR.l*((NC/vl - z)^2 * (1 - pl*LR.l) - NC/vl^2))
                                          # value with internal vl
      hess[ivl,ivl] <- hess[ivl,ivl]*vl^2 + grad[ivl]
                                          # Normalization for external vl
      hess[ipl,ipl] <- -sum(LR.l^2 - 2*LR.l*LR.h + LR.h^2)
      eb <- exp(-beta[ipl])
      dp2.db2 <- eb*(eb - 1)/(1 + eb)^3
      hess[ipl,ipl] <- hess[ipl,ipl]*(pl*(1 - pl))^2 + grad.ipl*dp2.db2
      hess[ipl,ivl] <- sum(LR.l*(NC/vl - z)*(1 - pl*(LR.l - LR.h)))
      hess[ipl,ivl] <- hess[ipl,ivl]*pl*(1 - pl)*vl
                                          # normalization
      hess[ivl,ipl] <- hess[ipl,ivl]
                                          # a-s are MxM matrices
      a1 <- t(dLl.dlL) %*% (NC*di - vl*mu*hazard.pieces)
      a2 <- -diag(as.vector(t(LR.l*vl*mu*hazard.pieces) %*% rep(1, N)))
      a3 <- t(dLh.dlL) %*% (NC*di - mu*hazard.pieces)
      a4 <- -diag(as.vector(t(LR.h*mu*hazard.pieces) %*% rep(1, N)))
                                          # and here we do a diagonal matrix from
                                          # the diagonal
      hess[ilambda,ilambda] <- pl*(a1 + a2) + ph*(a3 + a4) -
         t(vgrad.lambda) %*% vgrad.lambda
      if(!is.null(X)) {
         a1 <- pl*LR.l*((NC - z.l)^2 - z.l) + ph*LR.h*((NC - z)^2 - z) -
            vgrad.gamma^2
         hess[igamma,igamma] <- t(a1*X) %*% X
         a1 <- t(pl*dLl.dlL + ph*dLh.dlL) %*% (X*vgrad.gamma)
         a2 <- (NC - z.l)*(NC*di - vl*mu*hazard.pieces) - vl*mu*hazard.pieces
         a3 <- (NC - z)*(NC*di - mu*hazard.pieces) - mu*hazard.pieces
         a4 <- pl*t(a2) %*% (X*LR.l) + ph*t(a3) %*% (X*LR.h)
         hess[ilambda,igamma] <- a4 - a1
         hess[igamma,ilambda] <- t(hess[ilambda,igamma])
      }
      a1 <- pl*LR.l*(NC/vl - z)*vgrad.lambda
      a2 <- pl*LR.l*(NC/vl - z)*(NC*di - vl*mu*hazard.pieces)
      a3 <- pl*LR.l*mu*hazard.pieces
      hess[ilambda,ivl] <- t(-a1 + a2 - a3) %*% rep(1, N)
      hess[ilambda,ivl] <- hess[ilambda,ivl]*vl
                                          # normalization
      hess[ivl,ilambda] <- t(hess[ilambda,ivl])
      a1 <- (LR.l - LR.h)*vgrad.lambda
      a2 <- dLl.dlL - dLh.dlL
      hess[ilambda,ipl] <- t(-a1 + a2) %*% rep(1, N)
      hess[ilambda,ipl] <- hess[ilambda,ipl]*pl*(1 - pl)
                                          # normalization
      hess[ipl,ilambda] <- t(hess[ilambda,ipl])
      if(!is.null(X)) {
         a1 <- pl*LR.l*(NC/vl - z)*vgrad.gamma
         a2 <- pl*LR.l*(NC/vl - z)*(NC - z.l)
         a3 <- pl*LR.l*z
         hess[igamma,ivl] <- t(X) %*% (-a1 + a2 - a3)
         hess[igamma,ivl] <- hess[igamma,ivl]*vl
         hess[ivl,igamma] <- t(hess[igamma,ivl])
         a1 <- (LR.l - LR.h)*vgrad.gamma
         a2 <- LR.l*(NC - z.l) - LR.h*(NC - z)
         hess[igamma,ipl] <- t(X) %*% (-a1 + a2)
         hess[igamma,ipl] <- hess[igamma,ipl]*pl*(1 - pl)
                                          # normalization
         hess[ipl,igamma] <- t(hess[igamma,ipl])
      }
      ll <- -loglik
      attr(ll, "gradient") <- -grad
      attr(ll, "hessian") <- -hess
      ll
   }

   ##
   N <- length(t)
   if(cm[1] > 0) cm <- c(0, cm)
                                        # The periods start with time 0
   M <- length(cm)
   d <- ti(t, cm)
                                        # which hazard piece the time
                                        # falls to
   di <- matrix(0, N, M)
                                        # boolean indicator of exiting period
   di[cbind(1:N,d)] <- 1
   a <- matrix(0, M, M-1)
   a[lower.tri(a)] <- 1
   Di <- di %*% a
                                        # boolean for periods survived.  The last
                                        # period is never survived
   NC <- !censored
                                        # not censored
   if(!is.null(X)) {
      X <- as.matrix(X)
      K <- M + dim(X)[2] + 2
   } else {
      K <- M + 2
   }
                                        # number of parameters (2 for pl and vl)
   cat(N, "observations,", M, "time intervals and", K, "free parameters\n")
   cat(sum(NC), "non-censored and", sum(censored), "censored observations\n")
   ## now indexes of parameters
   ivl <- 1
   ipl <- ivl + 1
   ilambda <- (ipl + 1):(ipl + M)
   igamma <- (ipl + M + 1):K
   ##
   ## Now initial values
   if(is.null(b0)) {
      b0 <- c(-1, 0, rep(0,K-2))
      cm1 <- c(cm,+Inf)
                                        # here has last period end-value too
      l <- numeric(M)
      for(j in seq(cm)) {
         l[j] <- log(1/mean(t[t>cm1[j] & t<cm1[j+1]]))
                                        # initial value for lambda-terms is
                                        # average hazard in that time-piece
      }
      b0[ilambda] <- l
      cat("Initial value of parameter:\n")
      print(b0)
   }
   if(length(b0) != K) {
      stop(paste("length of parameter vector (currently ", length(b0), ")\n",
               "must equal the number of free parameters (",
               "currently ", K, ")", sep=""))
   }
   cat("Initial function value:", loglik(b0), "\n")
   a <- minNR(loglik, b0, ...)
   ## Now ensure that vl < 1.  If not, invert the coefficients and do just one
   ## iteration more
   if(normalise & a$code < 4 & a$estimate[ivl] > 1) {
      cat("Note: vl > 1.  Renormalizing to vl < 1\n")
      b0 <- a$estimate
      b0[ilambda] <- b0[ilambda] + b0[ivl]
      b0[ivl] <- -b0[ivl]
      b0[ipl] <- -b0[ipl]
      a1 <- minNR(loglik, b0, ...)
      a$maximum <- a1$maximum
      a$estimate <- a1$estimate
      a$gradient <- a1$gradient
      a$hessian <- a1$hessian
   }
   invisible(a)
}


