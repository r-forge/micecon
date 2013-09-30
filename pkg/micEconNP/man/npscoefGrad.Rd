\name{npscoefGrad}
\alias{npscoefGrad}

\title{Gradients for Smooth Coefficient Kernel Regression Models}

\description{
   Obtain gradients of the Z variables
   in non-parametric Smooth Coefficient Kernel Regression models
   by numeric finite-difference differenciation.
}

\usage{
   npscoefGrad( x, eps = 1e-3, \dots )
}

\arguments{
   \item{x}{an object of class \code{smoothcoefficient}
      (returned by \code{\link[np]{npscoef}}).}
   \item{eps}{the step size for the numeric finite-difference differentiation.}
   \item{\dots}{further arguments are currently ignored.}
}


\value{
   a three-dimensional array that contains the gradients
   of the (smooth) coefficients and the predicted values
   with respect to the Z variables.
   The first dimension of the array corresponds to the observations.
   The second dimension of the array corresponds to the (smooth) intercept,
   the (smooth) coefficients of the X variables, and
   the predicted values of the dependent (Y) variable
   (in that order).
   The third dimension of the array corresponds to the Z variables.
}

\seealso{\code{\link[np]{npscoef}}.}

\author{Arne Henningsen 
   (with significant help from Jeff Racine 
   who provided the idea of obtaining these gradients)}

\examples{
   n <- 250
   x1 <- runif(n)
   x2 <- rnorm(n)
   z1 <- runif(n)
   z2 <- runif(n, min=-2, max=2)
   y <- sin( 5 * z1 ) + x1 * exp( z2 ) + x2 * ( z1^2 + z2 ) + rnorm( n, sd = 0.2 )
   model <- npscoef( y ~ x1 + x2 | z1 + z2, beta = TRUE )
   
   print( model )
   summary( model )
   round( coef( model ), 3 )
   
   grad <- npscoefGrad( model )
   round( grad, 3 )
   
   all.equal( grad[ , 1, ] + grad[ , 2, ] * x1 + grad[ , 3, ] * x2, grad[ , 4, ] )
   
   \dontrun{
      for( j in 1:2 ) {
         for( i in 1:3 ) {
            plot( zData[,j], grad[ , i, j ],
               main = paste( "d b", i-1, " / d z", j, sep = "" ) )
            points(  zData[,j], gradTrue[ , i, j ], col = "red" )
            cat( "Press <enter>\n")
            readLines( n = 1 )
         }
      }
   }
}

\keyword{models}
\keyword{nonparametric}
\keyword{regression}
