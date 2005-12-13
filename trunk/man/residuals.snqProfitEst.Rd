\name{residuals.snqProfitEst}
\alias{residuals.snqProfitEst}

\title{Residuals of an SNQ profit function}

\description{
   Extract the residuals from the estimation of a Symmetric Normalized
   Quadratic (SNQ) profit function.
}

\usage{
\method{residuals}{snqProfitEst}( object, scaled = TRUE, ... )
}

\arguments{
   \item{object}{an object of type \code{snqProfitEst}.}
   \item{scaled}{logical. Return scaled quantities?}
   \item{...}{currently not used.}
}

\value{
   \code{residuals.snqProfitEst} returns a dataframe that
   contains the residuals for each netput and the profit.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}


\seealso{
   \code{\link{snqProfitEst}} and \code{\link{residuals}}
}

\examples{
   data( germanFarms )
   germanFarms$qOutput   <- germanFarms$vOutput / germanFarms$pOutput
   germanFarms$qVarInput <- -germanFarms$vVarInput / germanFarms$pVarInput
   germanFarms$qLabor    <- -germanFarms$qLabor
   germanFarms$time      <- c( 0:19 )
   pNames <- c( "pOutput", "pVarInput", "pLabor" )
   qNames <- c( "qOutput", "qVarInput", "qLabor" )
   estResult <- snqProfitEst( pNames, qNames, c("land","time"), data=germanFarms )
   residuals( estResult )
   residuals( estResult, scaled = FALSE )
}

\keyword{models}



