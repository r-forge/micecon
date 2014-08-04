library( "micEconNP" )
data( "oecdpanel" )

oecdpanel$yearFactor <- factor( oecdpanel$year )
npModel <- npreg( growth ~ yearFactor + initgdp, 
   regtype = "ll", gradients = TRUE, data = oecdpanel )

summary( npModel )

summary( npModel$bws )

bw <- sfactor2bw( npModel$bws$sfactor$x, c( "yearFactor", "initgdp" ),
   data = oecdpanel )
bw

all.equal( npModel$bws$bw, bw )
