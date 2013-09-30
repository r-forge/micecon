\name{npregGradFactor}
\alias{npregGradFactor}

\title{Gradients of Factor Variables in Kernel Regression Models}

\description{
   Obtain gradients of a categorical (\sQuote{factor}) variable
   in nonparametric kernel regression models.
}

\usage{
   npregGradFactor( x, varName )
}

\arguments{
   \item{x}{an object of class \code{npregression}
      (returned by \code{\link[np]{npreg}}).}
   \item{varName}{name of the categorical variable.}
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
   colMeans( yearGrad, na.rm = TRUE )
   colMedians( yearGrad, na.rm = TRUE )
}

\keyword{models}
\keyword{nonparametric}
\keyword{regression}
