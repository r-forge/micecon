## return the sum of an array while keeping its attributes
sumKeepAttr <- function( x, na.rm = FALSE ) {
   xAttr <- attributes( x )
   x <- sum( x, na.rm = na.rm )
   mostattributes( x ) <- xAttr
   return( x )
}
