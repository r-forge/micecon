\name{plotCrs}
\alias{plotCrs}

\title{Plot Results of CRS Models}

\description{
   This command creates plots
   that show the relationship between each explanatory variable
   and the dependent variable of CRS models
   holding all other explanatory variables
   constant at their median/modal values.
}

\usage{
   plotCrs( x )
}

\arguments{
   \item{x}{an object of class \code{crs}
      (returned by \code{\link[crs]{crs}}).}
}


\seealso{\code{\link[crs]{crs}}.}

\author{Arne Henningsen}

\examples{
   set.seed(42)
   n <- 250
   num.eval <- 50
   myData <- data.frame( x1 = runif( n ), x2 = runif( n ) )
   myData$z <- round( runif( n, min = 0, max = 3 ) )
   myData$dgp <- cos( 2 * pi * myData$x1 ) + sin( 2 * pi * myData$x2 ) +
      myData$z / 5
   myData$z <- factor( myData$z )
   myData$y <- myData$dgp + rnorm( n, sd = 0.5 )
   
   model <- crs( y ~ x1 + x2 + z, deriv = 1, nmulti = 1, data = myData )
   
   plotCrs( model )
}

\keyword{models}
\keyword{nonparametric}
\keyword{regression}
