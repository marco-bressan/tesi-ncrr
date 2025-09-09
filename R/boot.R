crr.boot <- function(data, stat, R, ran.gen, mle = NULL, ...,
                     parallel = FALSE, nclus = NA, trace = 1, seed = NULL) {

  para <- .setup.parallel(parallel, nclus, trace, seed, R,
                          pb = txtProgressBar(style = 3),
                          trace.init.msg = "Starting parametric bootstrap simulation")
  t0 <- stat(data, ...)
  res <- snowFT::performParallel(
    para$nclus,
    x = seq_len(R), initfun = \() devtools::load_all(), exitfun = \() message("Exit worker"),
    fun = function (i) {
      dd <- ran.gen(data, mle)
      stat(dd, ...)
    },
    printfun = para$printfun,
    printrepl = para$printrepl,
    ft_verbose = trace > 1,
    seed = para$seed
  )
  list(
    t0 = t0,
    t = do.call(rbind, res),
    R = R,
    data = data,
    statistic = stat,
    ran.gen = ran.gen,
    mle = mle,
    seed = para$seed
  )
}
