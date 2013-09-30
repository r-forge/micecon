library( "micEconNP" )
data( "oecdpanel" )

oecdpanel$yearFactor <- factor( oecdpanel$year )

npModel <- npreg( growth ~ yearFactor + initgdp, 
   regtype = "ll", gradients = TRUE, data = oecdpanel )

round( gradients( npModel ), 3 )

yearGrad <- npregGradFactor( npModel, "yearFactor" )

round( yearGrad, 3 )

round( colMeans( yearGrad, na.rm = TRUE ), 4 )

round( colMedians( yearGrad, na.rm = TRUE ), 4 )

all.equal( gradients( npModel )[ , 1 ], rowSums( yearGrad, na.rm = TRUE ) )

yearGradAll <- npregGradFactor( npModel, "yearFactor", onlyOwnLevels = FALSE )

round( yearGradAll, 3 )

round( colMeans( yearGradAll ), 4 )

round( colMedians( yearGradAll ), 4 )

all.equal( yearGradAll[ !is.na( yearGrad ) ], yearGrad[ !is.na( yearGrad ) ] )
