\name{MishraCES}
\docType{data}
\alias{MishraCES}

\title{Mishra's (2006) CES data}

\description{
   The \code{MishraCES} data set contains artificial production data.
   It has 50 observations (e.g. firms, sectors, or countries).
}

\usage{data(MishraCES)}

\format{
   This data frame contains the following columns:
   \describe{
      \item{No}{Firm number.}
      \item{Y}{Output quantity.}
      \item{X1}{Quantity of first input.}
      \item{X2}{Quantity of second input.}
      \item{X3}{Quantity of third input.}
      \item{X4}{Quantity of fouth input.}
   }
}

\source{
   Mishra, SK (2006):
   A Note on Numerical Estimation of Sato's Two-Level CES Production Function
   MPRA Working Paper No. 1019,
   \url{http://mpra.ub.uni-muenchen.de/1019/}.
}

\examples{
   # load the data set
   data( "MishraCES" )
   
   # show mean values of all variables
   colMeans( MishraCES )

   # re-calculate the endogenous variable (see Mishra 2006)
   # coefficients of the nested CES function with 4 inputs
   b <- c( "gamma" = 200 * 0.5^(1/0.6), "delta_1" = 0.6, "delta_2" = 0.3, 
      "delta" = 0.5, "rho_1" = 0.5, "rho_2" = -0.17, "rho" = 0.6 )
   MishraCES$Y2 <- cesCalc( xNames = c( "X1", "X2", "X3", "X4" ), 
      data = MishraCES, coef = b, nested = TRUE )
   all.equal( MishraCES$Y, MishraCES$Y2 )
}

\keyword{datasets}
