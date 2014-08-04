\name{sfactor2bw}
\alias{sfactor2bw}

\title{Convert scale factors to bandwidths}

\description{
   Convert scale factors (relative measures of bandwidths)
   to (absolute) bandwidths.
}

\usage{
   sfactor2bw( sfactor, varNames, data, order = 2, pre0.60 = FALSE )
}

\arguments{
   \item{sfactor}{a vector of scale factors of the explanatory variables.}
   \item{varNames}{a vector of the names of the explanatory variables.}
   \item{data}{data.frame with the that includes the explanatory variables.}
   \item{order}{order of the kernel for continuous explanatory variables.}
   \item{pre0.60}{logical.
      If \code{FALSE} (default),
      the conversion is done by the formula
      that is used in the \pkg{np} package version 0.60 or later.
      If \code{TRUE},
      the conversion is done by the formula
      that is used in the \pkg{np} package before version 0.60.}
}


\value{
   A vector of bandwidths.
}

\seealso{\code{\link[np]{npreg}}.}

\author{Arne Henningsen}

\examples{
   data( "oecdpanel" )
   
   oecdpanel$yearFactor <- factor( oecdpanel$year )
   npModel <- npreg( growth ~ yearFactor + initgdp, 
      regtype = "ll", gradients = TRUE, data = oecdpanel )
   
   all.equal( npModel$bws$bw,
      sfactor2bw( npModel$bws$sfactor$x, c( "yearFactor", "initgdp" ),
         data = oecdpanel ) )
}

\keyword{models}
\keyword{nonparametric}
\keyword{regression}
