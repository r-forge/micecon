sfactor2bw <- function( sfactor, varNames, data, order = 2, pre0.60 = FALSE ) {
   
   if( length( sfactor ) != length( varNames ) ) {
      stop( "argument 'sfactor' and 'varNames' must have the same length" )
   }
   
   nVar <- length( sfactor )
   isNumVar <- sapply( data[ , varNames ], is.numeric )
   nNumVar <- sum( isNumVar )
   
   sdNum <- apply( as.matrix( data[ , varNames[ isNumVar ] ] ), 2, sd )
   madNum <- apply( as.matrix( data[ , varNames[ isNumVar ] ] ), 2, mad )
   madNum[ madNum <= 0 ] <- Inf
   if( pre0.60 ) {
      madNum[ ] <- Inf
   }
   iqrNum <- apply( as.matrix( data[ , varNames[ isNumVar ] ] ), 2, IQR ) / 1.349
   iqrNum[ iqrNum <= 0 ] <- Inf
   
   spreadNum <- pmin( sdNum, madNum, iqrNum )
#    print( all.equal( spreadNum, np:::EssDee( data[ , varNames[ isNumVar ] ] ) ) )
   
#    print(rbind(sdNum,madNum,iqrNum,spreadNum))
   
   result <- rep( NA, nVar )
   
   result[ isNumVar ] <- sfactor[ isNumVar ] * spreadNum *
      nrow( data )^( -1 / ( 2 * order + nNumVar ) )
   
   result[ !isNumVar ] <- sfactor[ !isNumVar ] *
      nrow( data )^( -2 / ( 2 * order + nNumVar ) )

   result[ !isNumVar & result > 1 ] <- 1
   
   return( result )
}
