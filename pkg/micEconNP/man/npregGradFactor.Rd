\name{npregGradFactor}
\alias{npregGradFactor}

\title{Gradients of Factor Variables in Kernel Regression Models}

\description{
   Obtain gradients of an ordered or unordered categorical 
   (\sQuote{factor}) variable
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
   (except for the first level) on the dependent variable.
   In case of an unordered categorical variable,
   the current level is compared to the first level.
   In case of an ordered categorical variable,
   the current level is compared to the previous level.
   Each row of the matrix corresponds to an observation,
   while each column corresponds to a level.
}

\seealso{\code{\link[np]{npreg}}.}

\author{Arne Henningsen}

\examples{
   data( "oecdpanel" )

   ##############  unordered factor  ###################
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
   
   ##############  ordered factor  ###################
   oecdpanel$yearOrdered <- ordered( oecdpanel$year )
   npModelOrdered <- npreg( growth ~ yearOrdered + initgdp, 
      regtype = "ll", gradients = TRUE, data = oecdpanel )
   round( gradients( npModelOrdered ), 3 )

   yearGradOrdered <- npregGradFactor( npModelOrdered, "yearOrdered" )
   round( yearGradOrdered, 3 )
   colMeans( yearGradOrdered, na.rm = TRUE )
   colMedians( yearGrad, na.rm = TRUE )
   all.equal( gradients( npModelOrdered )[ oecdpanel$year != 1965, 1 ], 
      rowSums( yearGradOrdered, na.rm = TRUE )[ oecdpanel$year != 1965 ] )

   yearGradOrderedAll <- npregGradFactor( npModelOrdered, "yearOrdered", 
      onlyOwnLevels = FALSE )
   round( yearGradOrderedAll, 3 )
   colMeans( yearGradOrderedAll )
   colMedians( yearGradOrderedAll )

   all.equal( yearGradOrderedAll[ !is.na( yearGradOrdered ) ], 
      yearGradOrdered[ !is.na( yearGradOrdered ) ] )
}

\keyword{models}
\keyword{nonparametric}
\keyword{regression}
