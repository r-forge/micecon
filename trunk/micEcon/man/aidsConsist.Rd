\name{aidsConsist}
\alias{aidsConsist}
\alias{print.aidsConsist}
\alias{checkConsist.aidsEst}
\title{Check Consistency of the AIDS}
\description{
   Check whether the specification of the AIDS is consistent
   with microeconomic demand theory (i.e. utility maximisation).
}

\usage{
aidsConsist( priceNames, totExpName, coef, data,
   priceIndex = "TL", basePrices = NULL, baseShares = NULL,
   shareNames = NULL )

\method{print}{aidsConsist}( x, ... )

\method{checkConsist}{aidsEst}( object, ... )
}

\arguments{
   \item{priceNames}{a vector of strings containing the names of the prices.}
   \item{totExpName}{a string containing the variable name of total expenditure.}
   \item{coef}{a list containing the coefficients alpha, beta, gamma,
      and (only for the translog price index) alpha0.}
   \item{data}{a data frame containing the data.}
   \item{priceIndex}{a character string specifying the price index
      (see \code{\link{aidsPx}})
      or a numeric vector providing the log values of the price index
      (passed to \code{\link{aidsMono}}.}
   \item{basePrices}{a vector specifying the base prices for the
      Paasche, Laspeyres, and Tornqvist price index
      (passed to \code{\link{aidsMono}}.}
   \item{baseShares}{a vector specifying the base expenditure shares for the
      Laspeyres, simplified Laspeyres, and Tornqvist index
      (passed to \code{\link{aidsMono}}.}
   \item{shareNames}{a vector of strings containing the names of the expenditure
      shares (passed to \code{\link{aidsConcav}}.}
   \item{x}{an object of class \code{aidsConsist}.}
   \item{object}{an object of class \code{aidsEst}.}
   \item{...}{currently not used.}
}

\details{
   \code{checkConsist.aidsEst} is a wrapper function to \code{aidsConsist}
   that extracts all relevant arguments for \code{aidsConsist}
   from \code{object}.
}

\value{
   \code{aidsConsist} returns a list of class \code{aidsConsist}
   that contains following elements:
   \item{addingUp}{logical. Is the addinp-up condition fulfilled?}
   \item{homogeneity}{logical. Is the homogeneity condition fulfilled?}
   \item{symmetry}{logical. Is the symmetry condition fulfilled?}
   \item{mono}{an object returned by \code{\link{aidsMono}}.}
   \item{concav}{an object returned by \code{\link{aidsMono}}
      (only for the AIDS, i.e. if argument \code{priceIndex} is \code{"TL"}.}
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

\seealso{\code{\link{aidsEst}}}

\examples{
   data( Blanciforti86 )
   # Data on food consumption are available only for the first 32 years
   Blanciforti86 <- Blanciforti86[ 1:32, ]

   priceNames <- c( "pFood1", "pFood2", "pFood3", "pFood4" )
   shareNames <- c( "wFood1", "wFood2", "wFood3", "wFood4" )

   estResult <- aidsEst( priceNames, shareNames, "xFood",
      data = Blanciforti86, method = "IL" )
   aidsConsist( priceNames, "xFood", data = Blanciforti86,
      coef = estResult$coef )

   # the same can be obtained in an easier way
   checkConsist.aidsEst( estResult )
}

\keyword{models}
