# This file contains some functions used to transform the data.
eps <- .Machine$double.eps

# Logit base-10 transformation.
logit10 <- function(x)
  log10((x + eps)/(1 - x + eps))

# Project x onto the interval [a,b].
project.onto.interval <- function(x, a, b)
  pmin(b,pmax(a,x))

# Quantile normalization.
qt.random.tie <- function (x) {
  x.rank = rank(x,ties.method = "random")
  return(qqnorm(x.rank,plot.it = FALSE)$x)
}

# Returns the base-10 sigmoid of x. It is the inverse of logit10(x).
sigmoid10 <- function (x)
  1/(1 + 10^(-x))
