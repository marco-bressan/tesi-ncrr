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

crr.get.sigma <- function(design, ...) {
  dd <- design$design
  stopifnot("Baseline != 0 ancora da implementare!" = sapply(dd, \(d) d[1] == 0))
  Sigma <- blockdiag(lapply(dd, \(d) crr.vcov.baseline0(..., design = d)))
  Sigma
}

crr.get.mu <- function(design, ...) {
  dd <- design$design
  stopifnot("Baseline != 0 ancora da implementare!" = sapply(dd, \(d) d[1] == 0))
  mu <- do.call(c, lapply(dd, \(d) crr.mean.baseline0(..., design = d)))
  mu
}

crr.get.Gamma <- function(design, ...) {
  diag(design$gamma)
}

#' Parametri iniziali per l'ottimizzatore della verosimiglianza NCRR
#'
#' @export
#' @method getInitial ncrr.design
#' @param object Oggetto di tipo `ncrr.design`
#' @param data Non usato
#' @param ... Non usati
#' @param eps Tolleranza numerica per lo 0
#'
getInitial.ncrr.design <- function(object, data, ..., eps = 1e-10) {
  if (!missing(data)) .NotYetUsed("data", error = FALSE)
  nd <- length(unique(do.call(c, object$design)))
  structure(c(alpha = rep(0, nd - 1),
              beta = rep(1, nd - 1),
              mu0 = 0, sigma20 = 1, rho = 0,
              sigma = rep(1, nd - 1)),
            lower = c(alpha = rep(-Inf, nd - 1),
                      beta = rep(-Inf, nd - 1),
                      mu0 = -Inf, sigma20 = eps, rho = -1,
                      sigma = rep(eps, nd - 1)),
            upper = c(alpha = rep(Inf, nd - 1),
                      beta = rep(Inf, nd - 1),
                      mu0 = Inf, sigma20 = Inf, rho = 1,
                      sigma = rep(Inf, nd - 1)))
}
