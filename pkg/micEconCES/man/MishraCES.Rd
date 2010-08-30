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
}

\keyword{datasets}
