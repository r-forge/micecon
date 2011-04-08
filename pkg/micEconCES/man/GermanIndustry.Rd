\name{GermanIndustry}
\docType{data}
\alias{GermanIndustry}
\encoding{latin1}

\title{Aggregated Time Series Data for the West German Industry}

\description{
   The data frame \code{GermanIndustry} contains annual aggregated data
   of the entire West German industry from 1960 until 1993
   as well as data of seven industrial sectors from 1970 to 1988/1992.
   This data set has been used by Kemfert (1998).
}

\usage{data(GermanIndustry)}

\format{
   This data frame contains the following columns/variables:
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
      \item{C_Y}{gross value added of the West German
         chemical industry (in billion Deutsche Mark at prices of 1991).}
      \item{C_K}{capital: gross stock of fixed assets of the West German 
         chemical industry (in billion Deutsche Mark at prices of 1991).}
      \item{C_A}{labor: total number of persons employed in the West German 
         chemical industry (in thouands).}
      \item{C_E}{final energy consumption in the West German 
         chemical industry (in GWh).}
      \item{S_Y}{gross value added of the West German
         stone and earth industry (in billion Deutsche Mark at prices of 1991).}
      \item{S_K}{capital: gross stock of fixed assets of the West German 
         stone and earth industry (in billion Deutsche Mark at prices of 1991).}
      \item{S_A}{labor: total number of persons employed in the West German 
         stone and earth industry (in thouands).}
      \item{S_E}{final energy consumption in the West German 
         stone and earth industry (in GWh).}
      \item{I_Y}{gross value added of the West German
         iron industry (in billion Deutsche Mark at prices of 1991).}
      \item{I_K}{capital: gross stock of fixed assets of the West German 
         iron industry (in billion Deutsche Mark at prices of 1991).}
      \item{I_A}{labor: total number of persons employed in the West German 
         iron industry (in thouands).}
      \item{I_E}{final energy consumption in the West German 
         iron industry (in GWh).}
      \item{N_Y}{gross value added of the West German
         non-ferrous industry (in billion Deutsche Mark at prices of 1991).}
      \item{N_K}{capital: gross stock of fixed assets of the West German 
         non-ferrous industry (in billion Deutsche Mark at prices of 1991).}
      \item{N_A}{labor: total number of persons employed in the West German 
         non-ferrous industry (in thouands).}
      \item{N_E}{final energy consumption in the West German 
         non-ferrous industry (in GWh).}
      \item{V_Y}{gross value added of the West German
         vehicle industry (in billion Deutsche Mark at prices of 1991).}
      \item{V_K}{capital: gross stock of fixed assets of the West German 
         vehicle industry (in billion Deutsche Mark at prices of 1991).}
      \item{V_A}{labor: total number of persons employed in the West German 
         vehicle industry (in thouands).}
      \item{V_E}{final energy consumption in the West German 
         vehicle industry (in GWh).}
      \item{P_Y}{gross value added of the West German
         paper industry (in billion Deutsche Mark at prices of 1991).}
      \item{P_K}{capital: gross stock of fixed assets of the West German 
         paper industry (in billion Deutsche Mark at prices of 1991).}
      \item{P_A}{labor: total number of persons employed in the West German 
         paper industry (in thouands).}
      \item{P_E}{final energy consumption in the West German 
         paper industry (in GWh).}
      \item{F_Y}{gross value added of the West German
         food industry (in billion Deutsche Mark at prices of 1991).}
      \item{F_K}{capital: gross stock of fixed assets of the West German 
         food industry (in billion Deutsche Mark at prices of 1991).}
      \item{F_A}{labor: total number of persons employed in the West German 
         food industry (in thouands).}
      \item{F_E}{final energy consumption in the West German 
         food industry (in GWh).}
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
   Energy Economics 20: 249-264
   (doi:10.1016/S0140-9883(97)00014-5)
}

\keyword{datasets}

