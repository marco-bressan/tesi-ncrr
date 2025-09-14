crr.boot <- function(data, stat, R, ran.gen, mle = NULL, ...,
                     parallel = FALSE, nclus = NA, trace = 1, seed = NULL) {

  para <- .setup.parallel(parallel, nclus, trace, seed, R,
                          pb = txtProgressBar(style = 3),
                          trace.init.msg = "Starting parametric bootstrap simulation")
  t0 <- stat(data, ...)
  force(ran.gen)
  type0 <- typeof(t0)
  simple.type <- type0 %in% c("logical", "integer", "double", "complex", "character", "raw")
  len0 <- length(t0)
  res <- snowFT::performParallel(
    para$nclus,
    export = ls(envir = parent.frame()),
    x = seq_len(R), initfun = \() devtools::load_all(), exitfun = \() message("Exit worker"),
    fun = function(i) try({
      dd <- ran.gen(data, mle)
      stat(dd, ...)
    }),
    printfun = para$printfun,
    printrepl = para$printrepl,
    ft_verbose = trace > 1,
    seed = para$seed
  )
  close.parallel(para)
  if (simple.type) {
    errids <- which(sapply(res, \(x) inherits(x, "try-error")))
    if (length(errids) > 0) {
      conversion <- match.fun(paste0("as.", type0))
      for (i in errids) {
        tmp <- rep_len(cc <- conversion(res[[i]]), len0)
        attributes(tmp) <- attributes(res[[i]])
        if (length(cc) == 1 && is.na(cc))
          attr(tmp, ".error") <- paste(res[[i]], collapse = " ")
        res[[i]] <- tmp
      }
    }
  }

  out <- list(
    t0 = t0,
    t = if (simple.type) do.call(rbind, res) else res,
    R = R,
    data = data,
    statistic = stat,
    ran.gen = ran.gen,
    mle = mle,
    seed = para$seed
  )
  # copia degli attributi
  attr0 <- c(".error", names(attributes(t0)))
  for (i in setdiff(attr0, "names")) {
    aa <- lapply(res, \(x) attr(x, i))
    if (length(w <- which(sapply(aa, length) > 0)) > 0)
      attr(out, i) <- aa[w]
  }
  out
}
