## ---- calculation of netput quantities -------
snqProfitNetputs <- function( pNames, fNames, weights, coef ) {
   nNetput <- length( pNames )
   nFix    <- length( fNames )
   nObs <- nrow(P)
   X <- array(0,c(nObs,nNetput))
   w <- P%*%pw
   for(i in 1:nNetput) {
      X[,i] <- alpha[i]
      for(j in 1:nNetput) {
         X[,i] <- X[,i]+beta[i,j]*P[,j]/w
         for(k in 1:nNetput) X[,i] <- X[,i]-0.5*pw[i]*beta[j,k]*P[,j]*P[,k]/w^2
      }
      for(j in 1:nFix) {
         X[,i] <- X[,i]+delta[i,j]*Z[,j]
         for(k in 1:nFix) X[,i] <- X[,i]+0.5*pw[i]*gamma[j,k]*Z[,j]*Z[,k]
      }
   }
   # X <- array(1,c(nObs))%*%t(alpha)+(P%*%beta)/(w%*%array(1,c(1,nNetput)))-
   #         0.5*(diag(P%*%beta%*%t(P))%*%t(pw))/((w^2)%*%array(1,c(1,nNetput)))+
   #         Z%*%t(delta)+0.5*diag(Z%*%gamma%*%t(Z))%*%t(pw)
   return(X)
}
