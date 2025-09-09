.setup.parallel <- function(parallel, nclus, trace, seed, R, pb = NULL, trace.init.msg = NULL) {
  seed.in <- if (length(seed) == 0)
    sample.int(.Machine$integer.max, 6)
  else rep_len(as.numeric(seed), 6)
  # cluster setup
  printfun <- if (trace > 0)
    \(res, n, args = NULL) setTxtProgressBar(pb, n/length(res))
  if (isFALSE(parallel))
    nclus <- 0
  else if (anyNA(nclus))
    nclus <- parallel::detectCores()
  if (trace > 0) {
    cat(trace.init.msg,
        if(!isFALSE(parallel)) sprintf("(in parallel, %d core[s])", nclus),
        "\n")
    force(pb)
  }
  list(seed = seed.in, printfun = printfun, nclus = nclus,
       printrepl = max(min(50, R / 10), 1))
}

close.parallel <- function(con, ...) {
  if (is.null(con$printfun))
    invisible()
  env <- environment(con$printfun)
  if (exists("pb", envir = env) && !is.null(env$pb))
    close(env$pb)
  invisible()
}
