last <- function(v)
  v[length(v)]


tmapply <- function(X, INDEX, FUN=sum, ...) {
### tapply a function over a columns of a matrix
  tvapply <- function(X, ...) {
    ## apply the function for a (column) vector
    tapply(X, INDEX, FUN, ...)
  }
  apply(X, 2, FUN=tvapply, ...)
}
