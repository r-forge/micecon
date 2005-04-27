\name{invMillsRatio}
\alias{invMillsRatio}
\title{Calculates the 'Inverse Mills Ratio' of probit models}

\description{
   Calculates the 'Inverse Mills Ratio' of probit models.
}

\usage{
   invMillsRatio( x )
}

\arguments{
   \item{x}{probit model estimated by \code{\link{glm}} or
      \code{\link[VGAM]{vglm}}}
}

\value{
   A data frame that contains the Inverse Mills Ratio (IMR) and the delta
   values.

   If a univariate probit estimation is provided, the variables
   \code{IMR1} and \code{IMR0} are the Inverse Mills Ratios to correct
   for a sample selection bias of y = 1 and y = 0, respectively.
   Accordingly, 'delta1' and 'delta0' are the corresponding delta values.

   If a bivariate probit estimation is provided, the variables
   \code{IMRa1}, \code{IMRa0}, \code{IMRb1}, and \code{IMRb0} are the
   Inverse Mills Ratios to correct for a sample selection bias
   of y = 1 and y = 0 in equations 'a' and 'b', respectively.
   Accordingly, 'deltaa1', 'deltaa0', 'deltab1' and 'deltab0' are the
   corresponding delta values.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\examples{
## Wooldridge( 2003 ): example 17.5, page 590
data(Mroz87)
myProbit <- glm( lfp ~ nwifeinc + educ + exper + I( exper^2 ) + age +
   kids5 + kids618, family = binomial( link = "probit" ), data=Mroz87 )
Mroz87$IMR <- invMillsRatio( myProbit )$IMR1
myHeckit <- lm( log( wage ) ~ educ + exper + I( exper^2 ) + IMR,
   data = Mroz87[ Mroz87$lfp == 1, ] )

# using NO labor force participation as endogenous variable
Mroz87$nolfp <- 1 - Mroz87$lfp
myProbit2 <- glm( nolfp ~ nwifeinc + educ + exper + I( exper^2 ) + age +
   kids5 + kids618, family = binomial( link = "probit" ), data=Mroz87 )
all.equal( invMillsRatio( myProbit )$IMR1, invMillsRatio( myProbit2 )$IMR0 )
   # should be true
}

\keyword{models}
