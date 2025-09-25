.parallel <- function(x, fun, export = NULL, initfun = NULL, exitfun = NULL, ...,
                      parallel, nclus, trace, seed, pb = NULL,
                      trace.init.msg = NULL, min.repl = 10) {
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
  res <- snowFT::performParallel(nclus, export = export,
                                 x = x, initfun = initfun, exitfun = exitfun,
                                 fun = fun, printfun = printfun,
                                 printrepl = max(min(min.repl, length(x) / 10), 1),
                                 ft_verbose = trace > 1, seed = seed.in)
  attr(res, "seed") <- seed.in
  close(pb)
  res
}
