\name{utility}
\docType{data}
\alias{utility}

\title{Fossil Fuel Fired Steam Electric Power Plants}

\description{
   The data frame \code{utility} is a panel data set 
   of 72 investor-owned fossil fuel fired steam electric power plants
   in the United States over eleven years (1986-1996).
   This data set has 791 observations and, thus, is almost balanced
   (only one firm-year observation is missing).
   This data set is used as an example in 
   Kumbhakar, Wang, and Horncastle (2015).
}

\usage{data(utility)}

\format{
   This data frame contains the following columns/variables:
   \describe{
      \item{firm}{firm number (1-72).}
      \item{year}{last two digits of the year (86-96).}
      \item{y}{net steam electric power generation in MWh.}
      \item{regu}{a dummy variable indicating whether the firm is regulated.}
      \item{k}{capital input quantity.}
      \item{labor}{labor and maintenance input quantity.}
      \item{fuel}{fuel input quantity.}
      \item{wl}{price of the labour and maintenance input.}
      \item{wf}{price of the fuel input.}
      \item{wk}{price of the capital input.}
   }
}

\source{
   Companion website to Kumbhakar, Wang, and Horncastle (2015):
   \url{https://sites.google.com/site/sfbook2014/}.
}

\references{
   Kumbhakar, Subal C., Hung-Jen Wang, and Alan P. Horncastle (2015): 
   A Practitioner's Guide to Stochastic Frontier Analysis Using Stata. 
   Cambridge University Press. 
   (\href{http://doi.org/10.1017/CBO9781139342070}{doi:10.1017/CBO9781139342070})
}

\keyword{datasets}
