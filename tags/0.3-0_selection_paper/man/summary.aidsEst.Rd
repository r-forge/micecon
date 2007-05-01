\name{summary.aidsEst}
\alias{summary.aidsEst}
\alias{print.summary.aidsEst}
\title{Summarizing the Estimation of an Almost Ideal Demand System}

\description{
   \code{summary.aidsEst} summarizes the estimation results
   of an Almost Ideal Demand System (AIDS).
   Currently, it calculates demand elasticities at the sample mean values.
}

\usage{
\method{summary}{aidsEst}(object, elaFormula = NULL, quantNames = NULL, ... )

\method{print}{summary.aidsEst}( x, ... )
}

\arguments{
   \item{object}{an object of class \code{aidsEst}.}
   \item{elaFormula}{the formula to calculate the elasticities
      (see \code{\link{aidsEla}})
      If the demand system was estimated by the linear approximation (LA),
      the default formula is 'Ch'.
      If the demand system was estimated by the iterative linear least squares
      estimator (ILLE),
      the default formula is 'AIDS'.}
   \item{quantNames}{an optional vector of strings containing the names of
      the quantities (just to label elasticities).}
   \item{x}{an object of class \code{summary.aidsEst}.}
   \item{...}{currently ignored.}
}

\value{
   a list of class \code{summary.aidsEst} that contains
   all elements of the provided \code{object}
   plus following objects:
   \item{ela}{a list containing the elasticities
      (see \code{\link{aidsEla}}).}
   \item{elaFormula}{formula used to calculate the elasticities.}
}

\references{
   Deaton, A.S. and J. Muellbauer (1980)
   An Almost Ideal Demand System.
   \emph{American Economic Review}, 70, p. 312-326.

   Blundell, R. and J.M. Robin (1999)
   Estimationin Large and Disaggregated Demand Systems:
   An Estimator for Conditionally Linear Systems.
   \emph{Journal of Applied Econometrics}, 14, p. 209-232.
}

\seealso{\code{\link{aidsEst}}, \code{\link{aidsEla}}.}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\examples{
   # Using data published in Blanciforti, Green & King (1986)
   data( Blanciforti86 )
   # Data on food consumption are available only for the first 32 years
   Blanciforti86 <- Blanciforti86[ 1:32, ]

   ## Repeating the demand analysis of Blanciforti, Green & King (1986)
   estResult <- aidsEst( c( "pFood1", "pFood2", "pFood3", "pFood4" ),
      c( "wFood1", "wFood2", "wFood3", "wFood4" ), "xFood",
      data = Blanciforti86, method = "LA:SL",
      maxiter = 1, rcovformula = 1, tol = 1e-7 )
   summary( estResult, elaFormula = "Ch" )

   ## Repeating the evaluation of different elasticity formulas of
   ## Green & Alston (1990)
   priceNames <- c( "pFood1", "pFood2", "pFood3", "pFood4" )
   shareNames <- c( "wFood1", "wFood2", "wFood3", "wFood4" )

   # AIDS
   estResultA <- aidsEst( priceNames, shareNames, "xFood",
      data = Blanciforti86[ -1, ], maxiter = 1,
      rcovformula=1, tol=1e-7,
      method = "IL:L" )
   diag( summary( estResultA, elaFormula = "AIDS" )$ela$marshall )

   # LA-AIDS + formula of AIDS
   estResultL1 <- aidsEst( priceNames, shareNames, "xFood",
      data = Blanciforti86, maxiter = 1,
      rcovformula=1, tol=1e-7, method = "LA:SL" )
   diag( summary( estResultL1, elaFormula = "AIDS" )$ela$marshall )

   # LA-AIDS + formula of Eales + Unnevehr
   estResultL2 <- aidsEst( priceNames, shareNames, "xFood",
      data = Blanciforti86, maxiter = 1,
      rcovformula=1, tol=1e-7, method = "LA:SL" )
   diag( summary( estResultL2, elaFormula = "EU" )$ela$marshall )

   # LA-AIDS + formula of Chalfant:
   estResultL3 <- aidsEst( priceNames, shareNames, "xFood",
      data = Blanciforti86, maxiter = 1,
      rcovformula=1, tol=1e-7, method = "LA:SL" )
   diag( summary( estResultL3, elaFormula = "Ch" )$ela$marshall )

   # Estimations with a demand shifter: linear trend
   # LA-AIDS 
   Blanciforti86$trend <- c( 0:( nrow( Blanciforti86 ) - 1 ) )
   estResult <- aidsEst( priceNames, shareNames, "xFood",
      data = Blanciforti86, shifterNames = "trend" )
   summary( estResult )

   # AIDS
   estResult <- aidsEst( priceNames, shareNames, "xFood",
      data = Blanciforti86, shifterNames = "trend",
      maxiter = 1, method = "IL:L" )
   summary( estResult )

   # Estimations with two demand shifters: linear + quadratic trend
   # LA-AIDS
   Blanciforti86$trend2 <- c( 0:( nrow( Blanciforti86 ) - 1 ) )^2
   estResult <- aidsEst( priceNames, shareNames, "xFood",
      data = Blanciforti86, shifterNames = c( "trend", "trend2" ) )
   summary( estResult )

   # AIDS
   estResult <- aidsEst( priceNames, shareNames, "xFood",
      data = Blanciforti86, shifterNames = c( "trend", "trend2" ),
      maxiter = 1, method = "IL:L" )
   summary( estResult )
}

\keyword{models}
