#' Specifica il disegno della meta-analisi
#'
#'
#' @return lista con i dettagli dell'esperimento
#'
#' @param x un dataframe. Se contenente le colonne $m_{ik}$ e $s_{ik}$, queste
#'   vengono interpretate come medie e dev.std dello studio, che viene quindi
#'   assunto ad outcome normale. In questo caso,
#'   $\hat\Gamma=\text{diag}(s_{ik}^2)$ `x$treatment` dev'essere un fattore: il
#'   primo livello sarà il _baseline_.
#'
#' @export
ncrr.design <- function(x) {
  treatments <- levels(x$treatment)
  if (all(c("mik", "sik") %in% colnames(x))) {# risposta normale
    theta <- x$mik
    gamma <- x$sik^2
  } else {
    theta <- with(x, log(rik) - log(nik - rik))
    gamma <- with(x, 1 / rik + 1 / (nik + rik))
  }
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
  #stopifnot("Baseline != 0 ancora da implementare!" = sapply(dd, \(d) d[1] == 0))
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
  #stopifnot("Baseline != 0 ancora da implementare!" = sapply(dd, \(d) d[1] == 0))
  dunique <- unique(dd)
  #browser()
  mul <- lapply(dunique, \(d) {
    #if (!0 %in% d) browser()
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

##' L'uso della presente è confinato alle funzioni che calcolano i
##' parametri esplicitamente, contenute nella dir `data-raw/`
##'
##' @title Ottieni matrici dal design
##' @param object design di uno studio ncrr
##' @param what al momento solo "theta" o "gamma"
##' @return matrice di interesse (quella dei gamma o dei theta) in un
##'   formato che sia compatibile con la notazione usata in maxima nel calcolo
##'   delle derivate esplicite.
##' @author Marco Bressan
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
##' @param object Oggetto di tipo `ncrr.design`
##' @param data Non usato
##' @param ... Non usati
##' @param eps Tolleranza numerica per lo 0
##' @param seed Seme per la generazione di un punto di partenza casuale. Se si
##'   passa un singolo NA, il comportamento è completamente casuale e non
##'   riproducibile. Può essere un vettore con nomi dei parametri `c(alpha,
##'   beta, mu0, sigma20, rho, sigma2)`: per specificare quali NON
##'   inizializzare casualmente, mettere NA.
##' @param rep numero di random start
##' @param transform applica le trasformazioni ai parametri?
##' @param fixed elenco di nomi di parametri fissati (cioè esclusi dal vettore/matrice)
##' @param vcov.type preimpostazioni sui parametri fixed
##'
##' @return un vettore o una matrice con i punti di partenza casuali (se rep>1)
getInitial.ncrr.design <- function(object, data, ..., eps = 1e-10, seed = NA,
                                   rep = 1, transform = TRUE,
                                   fixed = NULL,
                                   vcov.type = c("normal", "simplified", "achana")) {
  vcov.type <- match.arg(vcov.type)
  if (!missing(data)) .NotYetUsed("data", error = FALSE)
  dd <- object$design
  nd <- length(unique(do.call(c, dd)))
  fixed <- union(fixed, switch(vcov.type, simplified = "sigma2",
                               achana = c("sigma2", "rho"),
                               default = NULL))
  parls <- crr.par.idx(nd - 1, fixed, lengths = TRUE)
  init <- mapply(rep, list(alpha = 0,
                           beta = 1,
                           mu0 = 0, sigma20 = 1 - isTRUE(transform), rho = 0,
                           sigma2 = 1 - isTRUE(transform))[names(parls)],
                 parls,
                 SIMPLIFY = FALSE)
  #browser()
  if (rep > 1) {
    #TODO: al momento rho è fuori scala. per semplicità si può operare sul logit e log per la varianza
    if (!is.null(names(seed))) {
      seed.old <- seed
      seed <- rep(NA_real_, length(init))
      names(seed) <- names(init)
      seed[names(seed.old)] <- seed.old
    } else if (length(seed) == 1) {
      if (!is.na(seed)) set.seed(seed)
      seed <- rpois(length(init), 400)
      names(seed) <- names(init)
    }
    noNAseed <- which(!is.na(seed))
    iidx <- match(names(noNAseed), names(init))
    #TODO: aggiungere switch per parametri univariati tipo mu e rho
    init[iidx] <- mapply(
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
      names(init)[iidx],
      grepl("sigma", names(init)[iidx], fixed = TRUE),
      SIMPLIFY = FALSE
    )

    #browser()
    init <- drop(do.call(cbind, init))
  } else {
    init <- do.call(c, init)
  }
  lower <- mapply(rep, list(alpha = -Inf, beta = -Inf,
                            mu0 = -Inf, sigma20 = eps, rho = -1,
                            sigma2 = eps)[names(parls)],
                  parls)
  upper <- mapply(rep, list(alpha = Inf, beta = Inf,
                            mu0 = Inf, sigma20 = Inf, rho = 1,
                            sigma2 = Inf)[names(parls)],
                  parls)
  if (!is.null(fixed)) {
    lower <- lower[-match(names(fixed), names(lower))]
    upper <- upper[-match(names(fixed), names(upper))]
  }
  structure(init,
            lower = unlist(lower),
            upper = unlist(upper))
}
