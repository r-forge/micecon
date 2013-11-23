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

##############  ordered factor  ###################
oecdpanel$yearOrdered <- ordered( oecdpanel$year )

npModelOrdered <- npreg( growth ~ yearOrdered + initgdp, 
   regtype = "ll", gradients = TRUE, data = oecdpanel )

round( gradients( npModelOrdered ), 3 )

yearGradOrdered <- npregGradFactor( npModelOrdered, "yearOrdered" )

round( yearGradOrdered, 3 )

round( colMeans( yearGradOrdered, na.rm = TRUE ), 4 )

round( colMedians( yearGrad, na.rm = TRUE ), 4 )

all.equal( gradients( npModelOrdered )[ oecdpanel$year != 1965, 1 ], 
   rowSums( yearGradOrdered, na.rm = TRUE )[ oecdpanel$year != 1965 ] )

yearGradOrderedAll <- npregGradFactor( npModelOrdered, "yearOrdered", 
   onlyOwnLevels = FALSE )

round( yearGradOrderedAll, 3 )

round( colMeans( yearGradOrderedAll ), 4 )

round( colMedians( yearGradOrderedAll ), 4 )

all.equal( yearGradOrderedAll[ !is.na( yearGradOrdered ) ], 
   yearGradOrdered[ !is.na( yearGradOrdered ) ] )


# test npregCv()
cv <- npregCv( npModel )
round( cv, 6 )
all.equal( cv, npModel$bws$fval )

cvOrdered <- npregCv( npModelOrdered )
round( cvOrdered, 6 )
all.equal( cvOrdered, npModelOrdered$bws$fval )
