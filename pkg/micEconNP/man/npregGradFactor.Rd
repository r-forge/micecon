\name{npregGradFactor}
\alias{npregGradFactor}

\title{Gradients of Factor Variables in Kernel Regression Models}

\description{
   Obtain gradients of a categorical (\sQuote{factor}) variable
   in nonparametric kernel regression models.
}

\usage{
   npregGradFactor( x, varName, onlyOwnLevels = TRUE )
}

\arguments{
   \item{x}{an object of class \code{npregression}
      (returned by \code{\link[np]{npreg}}).}
   \item{varName}{name of the categorical variable.}
   \item{onlyOwnLevels}{logical. If \code{TRUE} (default),
      the gradients with respect to a specific level
      are only retured for observations
      that are of this level.
      All other gradients are \code{NA}.}
}


\value{
   A matrix that contain the effects of each level of the chosen factor variable
   (except for the first level)
   on the dependent variable compared to the first level.
   Each row of the matrix corresponds to an observation,
   while each column corresponds to a level.
}

\seealso{\code{\link[np]{npreg}}.}

\author{Arne Henningsen}

\examples{
   data( "oecdpanel" )
   oecdpanel$yearFactor <- factor( oecdpanel$year )
   npModel <- npreg( growth ~ yearFactor + initgdp, 
      regtype = "ll", gradients = TRUE, data = oecdpanel )
   round( gradients( npModel ), 3 )
   yearGrad <- npregGradFactor( npModel, "yearFactor" )
   round( yearGrad, 3 )
   colMeans( yearGrad, na.rm = TRUE )
   colMedians( yearGrad, na.rm = TRUE )
   all.equal( gradients( npModel )[ , 1 ], rowSums( yearGrad, na.rm = TRUE ) )

   yearGradAll <- npregGradFactor( npModel, "yearFactor", onlyOwnLevels = FALSE )
   round( yearGradAll, 3 )
   colMeans( yearGradAll )
   colMedians( yearGradAll )

   all.equal( yearGradAll[ !is.na( yearGrad ) ], yearGrad[ !is.na( yearGrad ) ] )
}

\keyword{models}
\keyword{nonparametric}
\keyword{regression}
