library( micEcon )

## Missong77
data( Missong77 )
## price indices for Missong77
cat( "\nLaspeyres Price Indices for Missong77\n" )
print( priceIndex( c( "p.beef", "p.veal", "p.pork" ),
   c( "q.beef", "q.veal", "q.pork" ), 1, Missong77 ) )

cat( "\nPaasche Price Indices for Missong77\n" )
print( priceIndex( c( "p.beef", "p.veal", "p.pork" ),
   c( "q.beef", "q.veal", "q.pork" ), 1, Missong77, "Paasche" ) )

cat( "\nFisher Price Indices for Missong77\n" )
print( priceIndex( c( "p.beef", "p.veal", "p.pork" ),
   c( "q.beef", "q.veal", "q.pork" ), 1, Missong77, "Fisher" ) )

## quantity indices for Missong77
cat( "\nLaspeyres Quantity Indices for Missong77\n" )
print( quantityIndex( c( "p.beef", "p.veal", "p.pork" ),
   c( "q.beef", "q.veal", "q.pork" ), 1, Missong77 ) )

cat( "\nPaasche Quantity Indices for Missong77\n" )
print( quantityIndex( c( "p.beef", "p.veal", "p.pork" ),
   c( "q.beef", "q.veal", "q.pork" ), 1, Missong77, "Paasche" ) )

cat( "\nFisher Quantity Indices for Missong77\n" )
print( quantityIndex( c( "p.beef", "p.veal", "p.pork" ),
   c( "q.beef", "q.veal", "q.pork" ), 1, Missong77, "Fisher" ) )


## Bleymueller251
data( Bleymueller251 )
## price indices for Bleymueller251
cat( "\nLaspeyres Price Indices for Bleymueller251\n" )
print( priceIndex( c( "p.A", "p.B", "p.C", "p.D" ),
   c( "q.A", "q.B", "q.C", "q.D" ),  1, Bleymueller251 ) )

cat( "\nPaasche Price Indices for Bleymueller251\n" )
print( priceIndex( c( "p.A", "p.B", "p.C", "p.D" ),
   c( "q.A", "q.B", "q.C", "q.D" ), 1, Bleymueller251, "Paasche" ) )

cat( "\nFisher Price Indices for Bleymueller251\n" )
print( priceIndex( c( "p.A", "p.B", "p.C", "p.D" ),
   c("q.A", "q.B", "q.C", "q.D" ), 1, Bleymueller251, "Fisher" ) )

## quantity indices for Bleymueller251
cat( "\nLaspeyres Quantity Indices for Bleymueller251\n" )
print( quantityIndex( c( "p.A", "p.B", "p.C", "p.D" ),
   c("q.A", "q.B", "q.C", "q.D" ), 1, Bleymueller251 ) )

cat( "\nPaasche Quantity Indices for Bleymueller251\n" )
print( quantityIndex( c( "p.A", "p.B", "p.C", "p.D" ),
   c("q.A", "q.B", "q.C", "q.D" ), 1, Bleymueller251, "Paasche" ) )

cat( "\nFisher Quantity Indices for Bleymueller251\n" )
print( quantityIndex( c( "p.A", "p.B", "p.C", "p.D" ),
   c("q.A", "q.B", "q.C", "q.D" ), 1, Bleymueller251, "Fisher" ) )


## Blanciforti
data( Blanciforti86 )
## preparing data of Blanciforti
pNames <- c( "pFood1", "pFood2", "pFood3", "pFood4" )
qNames <- c( "qFood1", "qFood2", "qFood3", "qFood4" )
wNames <- c( "wFood1", "wFood2", "wFood3", "wFood4" )
for( i in 1:4 ) {
   Blanciforti86[[ qNames[ i ] ]] <- Blanciforti86[[ wNames[ i ] ]] *
      Blanciforti86[[ "xFood" ]] / Blanciforti86[[ pNames[ i ] ]]
}

## price indices for Blanciforti
cat( "\nLaspeyres Price Indices for Blanciforti\n" )
print( priceIndex( pNames, qNames, 1, Blanciforti86 ) )

cat( "\nPaasche Price Indices for Blanciforti\n" )
print( priceIndex( pNames, qNames, 1, Blanciforti86, "Paasche" ) )

cat( "\nFisher Price Indices for Blanciforti\n" )
print( priceIndex( pNames, qNames, 1, Blanciforti86, "Fisher" ) )

## quantity indices for Blanciforti
cat( "\nLaspeyres Quantity Indices for Blanciforti\n" )
print( quantityIndex( pNames, qNames, 1, Blanciforti86 ) )

cat( "\nPaasche Quantity Indices for Blanciforti\n" )
print( quantityIndex( pNames, qNames, 1, Blanciforti86, "Paasche" ) )

cat( "\nFisher Quantity Indices for Blanciforti\n" )
print( quantityIndex( pNames, qNames, 1, Blanciforti86, "Fisher" ) )


