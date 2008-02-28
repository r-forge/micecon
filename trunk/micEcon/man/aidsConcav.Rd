\name{aidsConcav}
\alias{aidsConcav}
\alias{print.aidsConcav}
\title{Concavity of the AIDS}
\description{
   Check whether the expenditure function of the AIDS is concave
   in prices.
}

\usage{
aidsConcav( priceNames, totExpName, coef, data, shareNames = NULL )

\method{print}{aidsConcav}( x, ... )
}

\arguments{
   \item{priceNames}{a vector of strings containing the names of the prices.}
   \item{totExpName}{a string containing the variable name of total expenditure.}
   \item{coef}{a list containing the coefficients alpha, beta, gamma,
      and (only for the translog price index) alpha0.}
   \item{data}{a data frame containing the data.}
   \item{shareNames}{a vector of strings containing the names of the expenditure
      shares.}
   \item{x}{an object of class \code{aidsConcav}.}
   \item{...}{currently unused.}
}

\details{
   If argument \code{shareNames} is specified,
   observed shares are used for the calculation of the 'C' matrices
   to check for concavity;
   if argument \code{shareNames} is \code{NULL} (e.g., not specified),
   fitted shares are used for the calculation of the 'C' matrices.
}

\value{
   \code{aidsConcav} returns a list of class \code{aidsConcav}
   that contains following elements:
   \item{concavity}{a logical vector indicating whether the concavity
      condition is fulfilled at each observation.}
   \item{cPercent}{percent of observations where the concavity condition is
      fulfilled.}
   \item{cMatrices}{a list of the 'C' matrices for each observation to check
      for concavity (see Deaton and Muellbauer, 1980b, p.76 ).}
}
\references{
   Deaton, A.S. and J. Muellbauer (1980a)
   An Almost Ideal Demand System.
   \emph{American Economic Review}, 70, p. 312-326.

   Deaton, A.S. and J. Muellbauer (1980b)
   \emph{Economics and Consumer Behavior},
   Cambridge University Press, Cambridge.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{\code{\link{aidsEst}}, \code{\link{aidsElas}}}

\examples{
   data( Blanciforti86 )
   # Data on food consumption are available only for the first 32 years
   Blanciforti86 <- Blanciforti86[ 1:32, ]

   priceNames <- c( "pFood1", "pFood2", "pFood3", "pFood4" )
   shareNames <- c( "wFood1", "wFood2", "wFood3", "wFood4" )

   estResult <- aidsEst( priceNames, shareNames, "xFood",
      data = Blanciforti86, method = "IL" )
   aidsConcav( priceNames, "xFood", coef = estResult$coef,
      data = Blanciforti86, shareNames = shareNames )
}

\keyword{models}
