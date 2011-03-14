\name{GermanIndustry}
\docType{data}
\alias{GermanIndustry}
\encoding{latin1}

\title{Aggregated Time Series Data for the West German Industry}

\description{
   The data frame \code{GermanIndustry} contains annual aggregated data
   of the West German industrial sector from 1960 until 1993.
   This data set has been used by Kemfert (1998).
}

\usage{data(GermanIndustry)}

\format{
   This data frame contains the following columns:
   \describe{
      \item{year}{the year.}
      \item{Y}{output: gross value added of the West German industrial sector
         (in billion Deutsche Mark at prices of 1991).}
      \item{K}{capital: gross stock of fixed assets of the West German 
         industrial sector (in billion Deutsche Mark at prices of 1991).}
      \item{A}{labor: total number of persons employed in the West German 
         industrial sector (in million).}
      \item{E}{energy: final energy consumption in the West German 
         industrial sector (in GWh).}
   }
}

\note{
   Please note that Kemfert (1998) disregards the years 1973-1975
   in her estimations due to economic disruptions.
}

\source{
   German Federal Statistical Office (Statistisches Bundesamt), 
   data taken from Kemfert (1998).
}

\references{
   Kemfert, Claudia (1998): 
   Estimated Substitution Elasticities of a Nested CES Production Funktion 
   Approach for Germany, 
   Energy Economics 20: 249-264.
}

\keyword{datasets}

