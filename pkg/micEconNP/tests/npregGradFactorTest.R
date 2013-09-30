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

