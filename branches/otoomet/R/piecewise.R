piecewise <- function(t, censored, X=NULL, cm, b0=NULL, ... ){
   ### t - measured time
   ### c - indicatiors for censoring 0 - uncensored, 1 - censored
   ### X - covariates w/o constant.  If NULL, then no covariates
   ### cm - start times for intervals (the first is 0 not to include)
   ### b0 - initial parameters (optional)
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

   loglik <- function(beta) {
      lambda <- beta[ilambda]
      gamma <- beta[igamma]
      if(!is.null(X)) {
         bX <- X %*% gamma
         mu <- exp(bX)
      } else {
         bX <- 0
         mu <- rep(1, N)
      }
      surv <- matrix(c(diff(cm),0), N, M, byrow=TRUE)*cbind(Di, 0)
      surv[cbind(1:N,d)] <- t - cm[d]
                                          # columns: time survived in particular
                                          # pieces.
      lambda.pieces <- matrix(exp(lambda), N, M, byrow=TRUE)
                                          # columns: corresponding baseline
                                          # hazard
      hazard.pieces <- surv*lambda.pieces
      z <- mu*((hazard.pieces) %*% rep(1,M))
                                          # survived time in the pieces *
                                          # corresponding baseline hazard, summed
                                          # over the pieces
      loglik <- sum(NC*(bX + lambda[d]) - z)
      ## now gradient
      grad <- numeric(K)
      grad[ilambda] <- t(NC*di) %*% rep(1,N) -
                                          # hazard part (if uncensored), sum
                                          # over all the individuals for all the
                                          # pieces
         t(hazard.pieces) %*% mu
                                          # time survived * corresponding hazard,
                                          # note that %*% sums over individuals
      if(!is.null(X)) {
         grad[igamma] <- t(NC - z) %*% X
      }
      ## now hessian
      hess <- matrix(NA, K, K)
      a <- -t(as.vector(mu)*hazard.pieces) %*% rep(1,N)
      hess[ilambda,ilambda] <- diag(as.vector(a))
      if(!is.null(X)) {
                                          # only if covariates exist
         hess[igamma,igamma] <- -t(as.vector(z)*X) %*% X
         hess[igamma,ilambda] <- -t(as.vector(mu)*X) %*% hazard.pieces
         hess[ilambda,igamma] <- t(hess[igamma,ilambda])
      }
      ll <- -loglik
      attr(ll, "gradient") <- -grad
      attr(ll, "hessian") <- -hess
      return( ll )
   }


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
      NX <- ncol(X)
   } else {
      NX <- 0
   }
   K <- M + NX
                                        # number of parameters
   cat(N, "observations,", M, "time intervals and", K, "free parameters\n")
   cat(sum(NC), "non-censored and", sum(censored), "censored observations\n")
   ## now indexes of parameters
   ilambda <- 1:M
   igamma <- (M+1):K

   if(is.null(b0)) {
      b0 <- c(rep(0,K))
      cm1 <- c(cm,+Inf)
                                        # here has last period end-value too
      l <- numeric(M)
      for(j in seq(cm)) {
         l[j] <- log(1/mean(t[t>cm1[j] & t<cm1[j+1]]))
                                        # initial value for lambda-terms is
                                        # average hazard in that time-piece
      }
      b0[ilambda] <- l
   }
   if(length(b0) != K) {
      stop(paste("length of parameter vector (currently ", length(b0), ")\n",
               "must equal the number of free parameters (",
               "currently ", K, ")", sep=""))
   }
   a <- minNR(loglik, b0, ...)
   invisible(a)
}
