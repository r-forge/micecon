\name{aidsElas}
\alias{aidsElas}
\alias{elas.aidsEst}
\alias{print.aidsElas}
\title{ Elasticities of the AIDS model }

\description{
   These functions calculate and print
   the demand elasticities of an AIDS model.
}

\usage{
   aidsElas( coef, shares, prices = NULL, formula = "AIDS",
      quantNames = NULL, priceNames = NULL, coefVcov = NULL, df = NULL )

   \method{elas}{aidsEst}( object, formula = NULL, ... )

   \method{print}{aidsElas}( x, ... )
}

\arguments{
   \item{coef}{a list containing the coefficients alpha, beta and gamma.}
   \item{shares}{a vector of the shares at which the elasticities should be calculated.}
   \item{prices}{a vector of the prices at which the elasticities should be calculated
      (only needed if formula 'AIDS' is used).}
   \item{formula}{the elsticity formula to be used (see details).}
   \item{quantNames}{an optional vector of strings containing the names of
      the quantities to label elasticities.}
   \item{priceNames}{an optional vector of strings containing the names of
      the prices to label elasticities.}
   \item{coefVcov}{variance covariance matrix of the coefficients (optional).}
   \item{df}{degrees of freedom to calculate P-values of the elasticities 
      (optional).}
   \item{object}{an object of class \code{aidsEst}.}
   \item{x}{an object of class \code{aidsElas}.}
   \item{...}{additional arguments of \code{elas.aidsEst}
      are passed to \code{aidsEla};
      additional arguments of \code{print.aidsElas}
      are currently ignored.}
}

\details{
   \code{elas.aidsEst} is a wrapper function to \code{aidsElas}
   that extracts the
   estimated coefficients (\code{coef}),
   mean expenditure shares (\code{wMeans}),
   mean prices (\code{pMeans}),
   names of the prices (\code{priceNames}),
   estimated coefficient variance covariance matrix (\code{coef$allcov}), and
   degrees of freedom (\code{est$df})
   from the object of class \code{aidsEst}
   and passes them to \code{aidsElas}.
   If argument \code{formula} in \code{elas.aidsEst} is not specified,
   the default value depends on the estimation method.
   If the demand system was estimated by the linear approximation (LA),
   the default formula is 'Ch'.
   If the demand system was estimated by the iterative linear least squares
   estimator (ILLE),
   the default formula is 'AIDS'.

   At the moment only the elasticity formulas of the orginal AIDS ('AIDS'),
   the formula of Chalfant (1987) ('Ch') and the formula of
   Eales and Unnevehr (1988) ('EU') are implemented.
   The variance covariance matrices of the elasticities are calculated using
   the formula of Klein (1953, p. 258) (also known as the delta method).
   At the moment this is implemented only for the elasticity formulas of the 
   orginal AIDS.
}

\value{
  a list of class \code{aidsElas} containing following elements:
  \item{formula}{the elasticity formula used to calculate these elasticities.}
  \item{exp}{vector of expenditure elasticities.}
  \item{hicks}{matrix of Hicksian (compensated) price elasticities.}
  \item{marshall}{matrix of Marshallian (uncompensated) price elasticities.}
  \item{allVcov}{variance covariance matrix of all elasticities.}
  \item{expVcov}{variance covariance matrix of the expenditure elasticities.}
  \item{hicksVcov}{variance covariance matrix of the Hicksian (compensated) 
    price elasticities.}
  \item{marshallVcov}{variance covariance matrix of the Marshallian 
    (uncompensated) price elasticities.}
  \item{expStEr}{standard errors of the expenditure elasticities.}
  \item{hicksStEr}{standard errors of the Hicksian (compensated) price 
    elasticities.}
  \item{marshallStEr}{standard errors of the Marshallian (uncompensated) 
    price elasticities.}
  \item{expTval}{t-values of the expenditure elasticities.}
  \item{hicksTval}{t-values of the Hicksian (compensated) price elasticities.}
  \item{marshallTval}{t-values of the Marshallian (uncompensated) price 
    elasticities.}
  \item{expPval}{P-values of the expenditure elasticities.}
  \item{hicksPval}{P-values of the Hicksian (compensated) price elasticities.}
  \item{marshallPval}{P-values of the Marshallian (uncompensated) price 
    elasticities.}
}

\references{
   Chalfant, J.A. (1987)
   A Globally Flexible, Almost Ideal Demand System.
   \emph{Journal of Business and Economic Statistics}, 5, p. 233-242.

   Deaton, A.S. and J. Muellbauer (1980)
   An Almost Ideal Demand System.
   \emph{American Economic Review}, 70, p. 312-326.

   Eales J.S. and L.J. Unnevehr (1988)
   Demand for beef and chicken products: separability and structural change.
   \emph{American Journal of Agricultural Economics}, 70, p. 521-532.
   
   Klein L.R. (1953)
   \emph{A Textbook of Econometrics}. Row, Petersen and Co., New York.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{\code{\link{aidsEst}}}

\examples{
   data( Blanciforti86 )
   # Data on food consumption are available only for the first 32 years
   Blanciforti86 <- Blanciforti86[ 1:32, ]

   estResult <- aidsEst( c( "pFood1", "pFood2", "pFood3", "pFood4" ),
      c( "wFood1", "wFood2", "wFood3", "wFood4" ), "xFood",
      data = Blanciforti86, method = "LA:L" )
   wMeans <- colMeans( Blanciforti86[ , c( "wFood1", "wFood2",
      "wFood3", "wFood4" ) ] )
   aidsElas( estResult$coef, wMeans, formula = "Ch" )

   ## Repeating the evaluation of different elasticity formulas of
   ## Green & Alston (1990)
   priceNames <- c( "pFood1", "pFood2", "pFood3", "pFood4" )
   shareNames <- c( "wFood1", "wFood2", "wFood3", "wFood4" )

   # AIDS
   estResultA <- aidsEst( priceNames, shareNames, "xFood",
      data = Blanciforti86[ -1, ],
      method = "IL:L" )
   diag( elas( estResultA, formula = "AIDS" )$marshall )

   # LA-AIDS + formula of AIDS
   estResultL1 <- aidsEst( priceNames, shareNames, "xFood",
      data = Blanciforti86, method = "LA:SL" )
   diag( elas( estResultL1, formula = "AIDS" )$marshall )

   # LA-AIDS + formula of Eales + Unnevehr
   estResultL2 <- aidsEst( priceNames, shareNames, "xFood",
      data = Blanciforti86, method = "LA:SL" )
   diag( elas( estResultL2, formula = "EU" )$marshall )

   # LA-AIDS + formula of Chalfant:
   estResultL3 <- aidsEst( priceNames, shareNames, "xFood",
      data = Blanciforti86, method = "LA:SL" )
   diag( elas( estResultL3, formula = "Ch" )$marshall )
}

\keyword{ models }
