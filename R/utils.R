blockdiag <- function(mats, fill = 0) {
  if (is.matrix(mats)) return(mats)
  stopifnot("only square matrices are supported!" = sapply(mats, nrow) == (nc <- sapply(mats, ncol)))
  out <- matrix(fill, sum(nc), sum(nc))
  ncc <- c(0, cumsum(nc))
  for (i in seq_along(mats)) {
    out[(ncc[i]+1):ncc[i + 1], (ncc[i]+1):ncc[i + 1]] <- mats[[i]]
  }
  out
}
