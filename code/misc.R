# This file contains functions that don't fit anywhere else.

# Creates a binary factor for each 'level' of the categorical variable
binary.from.categorical <- function (x, col.names) {
  y <- list()
  for (level in levels(x))
    y[[level]] <- factor(as.integer(x == level))
  out <- data.frame(y,check.names = FALSE)
  colnames(out) <- col.names
  return(out)
}
