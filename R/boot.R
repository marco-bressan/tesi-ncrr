crr.boot <- function(data, statistic, R, ran.gen, mle, retain.data = TRUE, ...,
                     parallel = FALSE, nclus = NA, trace = 1, seed = NULL) {
  opt.list <- list(...)
  R <- R[1] # per compatibilità con crr.boot.ci
  t0 <- statistic(data, ...)
  force(ran.gen)
  type0 <- typeof(t0)
  simple.type <- type0 %in% c("logical", "integer", "double", "complex", "character", "raw")
  len0 <- length(t0)
  boot.fun <- if (isTRUE(retain.data)) function(i) {
    dd <- ran.gen(data, mle)
    ss <- statistic(dd, ...)
    attr(ss, "data") <- dd
    attr(ss, "digest.data") <- digest::digest(dd)
    ss
  } else function(i) {
    dd <- ran.gen(data, mle)
    statistic(dd, ...)
  }
  res <- .parallel(
    parallel = parallel, nclus = nclus, trace = trace, seed = seed,
    pb = txtProgressBar(style = 3),
    trace.init.msg = "Starting parametric bootstrap simulation",
    export = ls(envir = parent.frame()),
    x = seq_len(R), initfun = \() devtools::load_all(),
    exitfun = \() message("Exit worker"),
    fun = boot.fun
  )
  seed <- attr(res, "seed")
  attr(res, "seed") <- NULL
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
    statistic = statistic,
    ran.gen = ran.gen,
    mle = mle,
    seed = seed,
    .dots = opt.list
  )
  class(out) <- "crr.boot"
  # copia degli attributi
  attr0 <- c(".error",
             if (isTRUE(retain.data)) "data",
             if (!isFALSE(retain.data)) names(attributes(t0)))
  for (i in setdiff(attr0, "names")) {
    aa <- lapply(res, \(x) attr(x, i))
    if (length(w <- which(sapply(aa, length) > 0)) > 0)
      attr(out, i) <- aa[w]
  }
  out
}

.poss <- function(t0, t, side, ...) {
  if (side == "lower")
    return(mean(t <= t0, ...))
  if (side == "upper")
    return(mean(t >= t0, ...))
  mean(t >= t0 | t <= t0, ...)
}

.mergelist <- function(l1, l2, nm = names(l2)) {
  for (p in nm) {
    if (p %in% names(l2))
      l1[[p]] <- l2[[p]]
  }
  l1
}

crr.boot.ci <- function(data, statistic, R, ran.gen, mle, ..., within = FALSE,
                        signif = 0.05, side = c("lower", "upper", "both"),
                        interval, grid.len,
                        psi.grid = seq(interval[1], interval[2], length.out = grid.len),
                        parallel = FALSE, nclus = NA, trace = 1, seed = NULL,
                        ignore.attrs = FALSE) {
  side <- match.arg(side)
  basic.args <- c("data", "statistic", "R", "ran.gen", "mle")
  dots <- boot.res <- list()
  if (inherits(data, "crr.boot")) {
    for (i in basic.args[-1])
      if (missingArg(as.symbol(i), eval = TRUE)) assign(i, data[[i]], environment())
    dots <- data$.dots
    if (!isTRUE(within)) {
      boot.res <- data
      boot.stat <- boot.res$t
      if (length(binf <- which(!is.finite(boot.stat))) > 0) {
        warning("Rimozione di ", length(binf), " valori della statistica non finiti.")
        boot.stat <- boot.stat[-binf]
      }
    }
    data <- data$data
  } else if (!isTRUE(within)) {
    boot.res <- crr.boot(data, statistic, R, ran.gen, mle, ...,
                         parallel = parallel, nclus = nclus, trace = trace,
                         seed = seed)
  }
  dots <- .mergelist(dots, list(...)) # i '...' sovrascrivono gli altri argomenti
  t0 <- {
    if (isTRUE(within))
      do.call(statistic, append(list(data), dots))
    else
      boot.res$t0
  }
  sarg2 <- names(formals(statistic))[1:2]
  xpars <- if (!isTRUE(ignore.attrs))
    setdiff(names(attributes(t0)), c("class", "names", "dim", "dimnames", ".error"))
  if (any((mm <- match(sarg2, names(dots), 0)) != 0 ))
    dots <- dots[-mm]
  arglist0 <- append(list(data, NA), attributes(t0)[xpars])
  r.val <- sapply(psi.grid, \(psi) {
    arglist0[[2]] <- psi
    do.call(statistic, append(arglist0, dots))
  })
  if (isFALSE(within)) {
    sigb2.val <- sapply(X = r.val, FUN = .poss, t = boot.res$t, side = side, na.rm = TRUE)
  } else {
    if (isTRUE(all.equal(seed, boot.res$seed))) {
      seed <- sample.int(.Machine$integer.max, 6)
      warning("Non è stato fornito esplicitamente il seed. ",
              "Impostato automaticamente a ", deparse1(seed))
    }
    # costruzione lista degli argomenti (NA è dove ci va psi[i])
    arglist <- arglist0
    # se presenti, si utilizzano gli stessi dati del campione bootstrap originale
    # e quindi la stessa stima mle in modo da velocizzare il processo.
    #TODO: si deve considerare come un metodo "quasi-esatto"?
    if (is.null(attr(boot.res, "data")))
      within <- TRUE
    dots2 <- dots
    if (any((mm <- match(names(dots), sarg2, 0)) != 0 ))
      dots2 <- dots[-mm]
    if (anyNA(within)) {
      arglist[[1]] <- attr(boot.res, "data")
      battrs <- attributes(boot.res)
      arglist <- .mergelist(arglist, battrs, xpars)
      ic.fun <- function(psi) {
        arglist[[2]] <- psi
        simplify2array(.mapply(statistic, arglist, MoreArgs = dots2))
      }
      ic.fun3 <- function(psi) {
        arglist[[2]] <- psi
        rr <- sample.int(1000, 2, replace = TRUE)
        sid <- \(n) (rr - 1) %% n + 1
        arglist2 <- lapply(arglist, \(x) x[sid(length(x))])
        r1 <- simplify2array(.mapply(statistic, arglist2, MoreArgs = dots2))
        arglist2 <- arglist2[1:2]
        r2 <- simplify2array(.mapply(statistic, arglist2, MoreArgs = dots2))
        if (!isTRUE(all.equal(r1, r2))) browser()
        r1
      }
      ic.fun2 <- function(psi) {
        arglist[[2]] <- psi
        res <- numeric(1000)
        for (i in seq_along(arglist[[1]])) {
          arglist.i <- lapply(arglist, \(x) x[[(i - 1) %% length(x) + 1]])
          ss <- do.call(statistic, append(arglist.i, dots2))
          if (anyNA(ss)) browser()
          res[i] <- ss
        }
        res
      }
    } else {
      names(arglist) <- c("data", sarg2[2])
      # argomenti per crr.boot
      dots2[["parallel"]] <- FALSE
      dots2[["retain.data"]] <- FALSE
      dots2[["statistic"]] <- statistic
      dots2[["R"]] <- R[if(length(R) > 1) 2 else 1]
      if (missing(ran.gen))
        stop("Non è stato fornito il generatore dei dati")
      dots2[["ran.gen"]] <- ran.gen
      dots2[["mle"]] <- mle
      dots2[["trace"]] <- 0
      ic.fun <- function(psi) {
        arglist[[2]] <- psi
        B <- do.call(crr.boot, append(arglist, dots2))
        c(B$t0, B$t)
      }
    }
    #browser()
    print(system.time(ic.fun(psi.grid[1])))
    inner.boot <- simplify2array(.parallel(
      parallel = parallel, nclus = nclus, trace = trace, seed = seed,
      pb = txtProgressBar(style = 3),
      trace.init.msg = "Calcolo intervalli bootstrap",
      export = ls(envir = parent.frame()),
      x = psi.grid, initfun = \() devtools::load_all(),
      exitfun = \() message("Exit worker"),
      fun = ic.fun
    ), higher = FALSE)
    save(inner.boot, file = "innerboot_exact.rda")
    if (anyNA(within))
      inner.boot <- rbind(r.val, inner.boot)
    sigb2.val <- apply(inner.boot, 2, \(x) .poss(x[1], x[-1], side, na.rm = TRUE))
  }
  browser()
  sigb2.val <- clamp(sigb2.val, eps = 1e-8)
  sm1 <- smooth.spline(qnorm(sigb2.val), psi.grid)
  # intervalli di confidenza
  ic.vals <- c(signif / 2, .5, 1 - signif / 2)
  if (side == "lower") ic.vals <- rev(ic.vals)
  ic <- predict(sm1, qnorm(ic.vals))[["y"]]
  structure(list(ic = ic, spline = sm1), class = "crr.boot.ci")
}
