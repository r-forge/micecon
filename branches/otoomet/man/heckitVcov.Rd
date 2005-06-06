\name{heckitVcov}
\alias{heckitVcov}
\title{Heckit Variance Covariance Matrix}

\description{
   Calculate Variance Covariance Matrix for Heckit Estimations.
}

\usage{
   heckitVcov( xMat, wMat, vcovProbit, rho, delta, sigma,
   saveMemory = TRUE )
}

\arguments{
   \item{xMat}{model matrix of the second step estimation.}
   \item{wMat}{model matrix of the first step probit estimation.}
   \item{vcovProbit}{variance covariance matrix of the first step 
      probit estimation.}
   \item{rho}{.}
   \item{delta}{.}
   \item{sigma}{.}
   \item{saveMemory}{logical. Save memory by using a different algorithm?
      (this should not influence the results)}
}

\value{
   the variance covariance matrix.
}

\seealso{\code{\link{heckit}}.}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\examples{
}

\keyword{models}
