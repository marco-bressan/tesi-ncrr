#' Specifica il disegno della meta-analisi
#'
#' **NOTA:** al momento bisogna assicurarsi che il dataframe sia ordinato per studio
#' e per trattamento, con il baseline in prima posizione
#'
#' @return lista con i dettagli dell'esperimento
#'
#' @param x un dataframe
#' @param baseline livello baseline (di default il primo)
#'
#' @export
ncrr.design <- function(x, baseline) {
  treatments <- levels(x$treatment)
  if (missing(baseline))
    baseline <- treatments[1]
  theta <- with(x, log(rik) - log(nik - rik))
  gamma <- with(x, 1 / rik + 1 / (nik + rik))
  design <- tapply(x$treatment, x$study.id, \(z) as.integer(z) - 1)
  if (any(sapply(design, \(d) d[1] != 0))) {
    warning("Alcuni studi non hanno baseline 0: le corrispondenti routine non sono al momento implementate!")
  }
  if (!is.character(x$study.id))
    names(design) <- NULL
  ll <- as.list(environment())
  class(ll) <- "ncrr.design"
  return(ll)
}


#' @export
#' @method print ncrr.design
print.ncrr.design <- function(x, ...) {
  cat("Una Network Meta-Analisi con ", length(x$design), "studi e", length(x$treatments), "trattamenti:\n\n")
  ff <- factor(x$treatments[do.call(c, x$design) + 1], levels = x$treatments)
  dd <- as.data.frame(table(ff))
  colnames(dd) <- c("Trattamento", "N. occorrenze negli studi")
  rownames(dd) <- names(x$design)

  print(dd)
  cat("\n")
}


#' @export
#' @method subset ncrr.design
subset.ncrr.design <- function(x, subset, ...) {
  if (is.null(idx <- names(x$design)))
    idx <- seq_along(x$design)
  subs <- stats::na.omit(match(subset, idx))
  if (anyNA(subs)) {
    warning("Forniti identificatori di studio sconosciuti: ",
            paste(unique(subset[attr(subs, "na.action")]), collapse = ", "))
  }
  idx <- which(rep(seq_along(x$design), lengths(x$design)) %in% subs)
  x$theta <- x$theta[idx]
  x$gamma <- x$gamma[idx]
  x$design <- x$design[subs]

  return(x)
}

crr.get.sigma <- function(object, ..., raw = FALSE,
                          type = c("normal", "simplified", "achana")) {
  dd <- object$design
  type <- match.arg(type)
  vcovfn <- switch(type,
                   "normal" = crr.vcov,
                   "simplified" = crr.vcov.simple,
                   "achana" = crr.vcov.achana)
  stopifnot("Baseline != 0 ancora da implementare!" = sapply(dd, \(d) d[1] == 0))
  dunique <- unique(dd)
  #browser()
  Sigmal <- lapply(dunique, \(d) {
    psel <- par.select.multi(d, ...)
    psel[["design"]] <- d
    do.call(vcovfn, psel)
  })
  if (raw)
    return(Sigmal[match(dd, dunique)])
  blockdiag(Sigmal[match(dd, dunique)])
}

crr.get.mu <- function(object, ..., raw = FALSE) {
  dd <- object$design
  stopifnot("Baseline != 0 ancora da implementare!" = sapply(dd, \(d) d[1] == 0))
  dunique <- unique(dd)
  #browser()
  mul <- lapply(dunique, \(d) {
    psel <- par.select.multi(d, ...)
    psel[["design"]] <- d
    do.call(crr.mean.baseline0, psel)
  })
  if (raw)
    return(mul[match(dd, dunique)])
  do.call(c, mul[match(dd, dunique)])
}

par.select.multi <- function(object, ...) {
  lapply(list(...), \(x) if (length(x) == 1) x else x[setdiff(object, 0)])
}

crr.get.Gamma <- function(object, ..., raw = FALSE) {
  if (!raw)
    return(diag(object$gamma))
  lc <- cumsum(ll <- lengths(object$design))
  lc <- c(0, lc[-length(lc)])
  mapply(\(pos, len) diag(object$gamma[(pos + 1):(pos + len)]),
         lc, ll, SIMPLIFY = FALSE)
}

crr.get.theta <- function(object, ..., raw = FALSE) {
  if (!raw)
    return(object$theta)
  lc <- cumsum(ll <- lengths(object$design))
  lc <- c(0, lc[-length(lc)])
  mapply(\(pos, len) object$theta[(pos + 1):(pos + len)],
         lc, ll, SIMPLIFY = FALSE)
}

get.matrix.from.design <- function(object, what = c("theta", "gamma")) {
  what <- match.arg(what, several.ok = FALSE)
  stopifnot("passato oggetto non valido" = length(ll <- lengths(object$design)) > 0,
            "alcuni design sono nulli" = ll != 0)
  maxd <- max(unlist(object$design))
  mm <- matrix(NA, maxd + 1, length(object$design))
  for (i in seq_along(ll)) {
    mm[object$design[[i]] + 1, i] <- 1
  }
  stopifnot("lunghezza di `what` non compatibile" =
              sum(mm, na.rm = TRUE) == length(object[[what]]))
  mm[!is.na(mm)] <- object[[what]]
  t(mm)
}


#' Parametri iniziali per l'ottimizzatore della verosimiglianza NCRR
#'
#' @export
#' @method getInitial ncrr.design
#' @param object Oggetto di tipo `ncrr.design`
#' @param data Non usato
#' @param ... Non usati
#' @param eps Tolleranza numerica per lo 0
#' @param seed Seme per la generazione di un punto di partenza casuale (se NA non
#' genera valori casuali). Può essere un vettore con alcuni o tutti i nomi
#' dei parametri `c(alpha, beta, mu0, sigma20, rho, sigma2)`
#' per specificare quali inizializzare casualmente e quali no.
#'
getInitial.ncrr.design <- function(object, data, ..., eps = 1e-10, seed = NA,
                                   rep = 1, transform = TRUE) {
  if (!missing(data)) .NotYetUsed("data", error = FALSE)
  dd <- object$design
  nd <- length(unique(do.call(c, dd)))
  init <- list(alpha = rep(0, nd - 1),
               beta = rep(1, nd - 1),
               mu0 = 0, sigma20 = 1 - isTRUE(transform), rho = 0,
               sigma2 = rep(1 - isTRUE(transform), nd - 1))
  #browser()
  if (any(!is.na(seed))){
    #TODO: al momento rho è fuori scala. per semplicità si può operare sul logit e log per la varianza
    if (!is.null(names(seed))) {
      seed.old <- seed
      seed <- rep(NA_real_, length(init))
      names(seed) <- names(init)
      seed[names(seed.old)] <- seed.old
    } else if (length(seed) == 1) {
      set.seed(seed)
      seed <- rpois(length(init), 400)
    }
    noNAseed <- which(!is.na(seed))
    #TODO: aggiungere switch per parametri univariati tipo mu e rho
    init[noNAseed] <- mapply(
      \(s, nm, is.var.term) {
        if (is.na(s)) return(rep(NA, length))
        if (length(seed) > 1) set.seed(s)
        dd <- rnorm(rep * (nd - 1), 0, 4)
        if (!isTRUE(transform) && is.var.term)
          dd <- dd^2/4
        if (nm %in% c("mu0", "sigma20", "rho"))
          return(dd[1:rep])
        M <- matrix(dd, rep, nd - 1)
        colnames(M) <- paste0(nm, if (nd > 2) seq_len(ncol(M)) else "")
        M
      },
      seed[noNAseed],
      names(init)[noNAseed],
      grepl("sigma", names(init)[noNAseed], fixed = TRUE),
      SIMPLIFY = FALSE
    )

    #browser()
    init <- drop(do.call(cbind, init))
  } else {
    init <- do.call(c, init)
  }
  structure(init,
            lower = if (!isTRUE(transform)) {
              c(alpha = rep(-Inf, nd - 1),
                beta = rep(-Inf, nd - 1),
                mu0 = -Inf, sigma20 = eps, rho = -1,
                sigma2 = rep(eps, nd - 1))
            },
            upper = if (!isTRUE(transform)) {
              c(alpha = rep(Inf, nd - 1),
                beta = rep(Inf, nd - 1),
                mu0 = Inf, sigma20 = Inf, rho = 1,
                sigma2 = rep(Inf, nd - 1))
            })
}
