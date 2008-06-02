\name{summary.front41Output}
\alias{summary.front41Output}
\title{Summarizing the Estimation of Frontier 4.1}

\description{
   \code{summary.front41Output} summarizes the estimation results
   of a model estimated by Frontier 4.1..
}

\usage{
\method{summary}{front41Output}( object, ... )
}

\arguments{
   \item{object}{an object of class \code{front41Output}
      (read/created by \code{\link{front41ReadOutput}}.}
   \item{...}{currently ignored.}
}

\value{
   The \code{summary} method returns a list of class
   \code{summary.front41Output}
   with the same elements as object returned by \code{\link{front41ReadOutput}}.
   However, the elements \code{olsResults}, \code{gridResults}, and
   \code{mleResults} have an additional culumn with marginal significance
   levels (P values).
   The \eqn{P} values of the OLS estimates are calculated using
   the \eqn{t} distribution,
   while the (asymptotic) \eqn{P} values of the ML estimates
   are calculated based on the assumption
   that their \eqn{t} values follow an (asymptotic) standard normal
   distribution.
}

\seealso{\code{\link{aidsEst}}, \code{\link{aidsElas}}.}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\examples{
   # Using data published in Blanciforti, Green & King (1986)
   data( Blanciforti86 )
   # Data on food consumption are available only for the first 32 years
   Blanciforti86 <- Blanciforti86[ 1:32, ]

   ## Repeating the demand analysis of Blanciforti, Green & King (1986)
   estResult <- aidsEst( c( "pFood1", "pFood2", "pFood3", "pFood4" ),
      c( "wFood1", "wFood2", "wFood3", "wFood4" ), "xFood",
      data = Blanciforti86, priceIndex = "SL" )
   print( summary( estResult ) )
}

\keyword{models}
