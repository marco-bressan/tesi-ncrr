blockdiag <- function(mats, fill = 0) {
  if (is.matrix(mats)) return(mats)
  nc <- sapply(mats, ncol)
  nr <- sapply(mats, nrow)
  out <- matrix(fill, sum(nr), sum(nc))
  ncc <- c(0, cumsum(nc))
  nrc <- c(0, cumsum(nr))
  for (i in seq_along(mats)) {
    out[(nrc[i]+1):nrc[i + 1], (ncc[i]+1):ncc[i + 1]] <- mats[[i]]
  }
  out
}
