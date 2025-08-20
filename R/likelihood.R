##' A partire da un *design*, crea un oggetto funzione che può essere
##' passato ad `optim`.
##'
##' @export
##' @title Funzione per la definizione della verosimiglianza per NCRR.
##' @param object oggetto che definisce il design
##' @param transform applica trasformazioni ai parametri di
##'   varianza/correlazione
##' @param echo regola il livello delle stampe di debug
##' @param vcov.type struttura della matrice di varianza-covarianza.
##' @return una funzione del vettore dei parametri (in tal senso `llik1` ne
##'   costituisce una versione semplificata), con la possibilità di fissare gli
##'   stessi (opzione `fixed = list(...)`)
##' @author Marco Bressan
get.llik.from.design <- function(object, transform = TRUE, echo = 0,
                                 vcov.type = attr(object, "vcov.type")) {
  np <- length(tt <- unique(do.call(c, object$design))) - 1
  fixed.default <- NULL
  if (is.null(vcov.type))
    vcov.type <- "normal"
  vcov.type <- match.vcov.type(vcov.type)
  fixed.default <- match.vcov.fixed(vcov.type, TRUE, np)
  GETPARS <- function(params, fixed) {
    if (length(names(fixed.default)) > 0) {
      fixed[names(fixed.default)] <- fixed.default
      attributes(fixed) <- attributes(fixed.default)
    }
    params <- crr.split.par(params, np, transform = transform,
                            fixed = names(fixed),
                            parlen = attr(fixed, "parlen"))
    if (length(fixed) > 0)
      params <- append(fixed, params)
    if (vcov.type != "normal")
      return(set.vcov.params(params, np, vcov.type))
    params
  }

  # questo è il return,  la funzione obiettivo
  function(params, y = crr.get.theta(object, raw = TRUE),
           Gamma = crr.get.Gamma(object, raw = TRUE),
           fixed = NULL) {
    params <- GETPARS(params, fixed)
    if (echo > 1) {
      mapply( \(x, nm) paste(nm, "=", deparse1(round(x, 6))),
             params, names(params)) |>
        paste(collapse = ", ") |>
        cat("\n")
    }
    mu <- crr.get.mu(object, params, raw = TRUE)
    Sigma <- crr.get.sigma(object, params, raw = TRUE)

    ll <- mapply(\(t, m, s, g) {
      mvtnorm::dmvnorm(t, mean = m, sigma = s + g, log = TRUE)
    }, y, mu, Sigma, Gamma)
    if (anyNA(ll)) {
      warning("Si sono prodotti NA nel calcolo della verosimiglianza,",
              " che sono stati scartati.")
      stop("rilevati NA nella verosimiglianza")
    }
    if (echo > 2) {
      cat("CURRENT PIECEWISE LLIK:\n")
      mapply(\(m, s, l) {
        colnames(s) <- c("SIGMA", rep("", ncol(s) - 1))
        print(cbind("MU" = m, s))
        print(c("LLIK" = l))
      }, mu, Sigma, ll)
    }
    #browser()
    ll <- sum(ll, na.rm = TRUE)
    if (echo > 0)
      cat("CURRENT VALUE: ", ll, "\n")
    if (echo > 1)
      cat("=====================================\n")
    return(ll)
  }
}
