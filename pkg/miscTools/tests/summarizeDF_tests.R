## testing summarizeDF()
library( "miscTools" )
data( "cars" )

summarizeDF( cars )

summarizeDF( cars, maxLevel = 10 )

summarizeDF( cars, printValues = FALSE )


data( "iris" )

summarizeDF( iris )

summarizeDF( iris, printValues = FALSE )


data( "swiss" )

summarizeDF( swiss )

summarizeDF( swiss, maxLevel = 10 )

summarizeDF( swiss, printValues = FALSE )
