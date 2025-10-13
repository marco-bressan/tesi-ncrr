crr.mu.int <- function(pp, baseline = TRUE) {
  if (baseline) return(c(0, pp$alpha) + c(1, pp$beta) * pp$mu0)
  # studi non baseline
  if (USA_MIA_MODELLAZIONE) {
    #message("Chiamata mia implementazione!")
    mu_ib <- pp$alpha[1] + pp$beta[1] * pp$mu0
    mu <- c(mu_ib, pp$alpha[-1] - pp$alpha[1] + (pp$beta[-1] - pp$beta[1]) * mu_ib)
    return(mu)
  }
  #browser()
  return(c(diff(pp$alpha),
           diff(rev(pp$alpha))) + c(diff(pp$beta), diff(rev(pp$beta))) * pp$mu0)
}

crr.sigma.int <- function(pp, baseline = TRUE) {
  if (baseline) {
    # minore matrice vcov ottenuto togliendo la prima riga e la prima colonna
    #if (length(beta) > 1 && all(c(beta, sigma20, sigma2) != 1)) browser()
    pp$sigma2 <- sqrt(pp$sigma2)
    #if (any(!is.finite(sigma2))) browser()#stop("sigma2 negativo!")
    vv <- tcrossprod(c(1, pp$beta)) * pp$sigma20
    vv[-1, -1] <- vv[-1, -1] + tcrossprod(pp$sigma2) *
      diagoffdiag(1, pp$rho, length(pp$beta))
    return(vv)
  }
  sigmab <- pp$rho * sqrt(pp$sigma2[-1] * pp$sigma2[1]) # c(sigma^2_12, sigma^2_21)
  if (USA_MIA_MODELLAZIONE) {
    # in questa parte del codice voglio provare ad implementare la mia
    # versione della ncrr senza baseline. si può cambiare settando la
    # variabile globale (a livello di pacchetto) pari a FALSE: in questo caso
    # si userà la parametrizzazione di Guolo
    # !!! SI ASSUME CHE il baseline SIA IN PRIMA POSIZIONE !!!
    # === DA RIVEDERE ===
    betab <- pp$beta[-1] - pp$beta[1]
    vub <- pp$beta[1]^2 * pp$sigma20 + pp$sigma2[1]
    vv <- tcrossprod(c(1, betab)) * vub
    vv[-1, -1] <- vv[-1, -1] +
      diagoffdiag(1, pp$rho, length(pp$beta) - 1) * # len(beta) == 2 per design bivar.
      tcrossprod(sigmab)
    return(vv)
  }
  beta <- c(pp$beta[1] - pp$beta[2], pp$beta[2] - pp$beta[1])
  # sigma_12 = cov(eps1_01, eps1_02) = rho * sigma_01 * sigma_02
  return(tcrossprod(beta) * pp$sigma20 + diag(sigmab, 2))
}

##' Funzioni a basso livello per il calcolo dei parametri della distribuzione
##' marginale (normale) della NCRR
##'
##' NOTA: le funzioni `crr.vcov.*()` sono state deprecate e non dovrebbero più
##' essere utilizzate
##'
##' @title Parametri per la verosimiglianza normale
##' @param params Il vettore dei parametri.
##' @param design Un vettore intero indicante il design dello studio in esame.
##' @return Un vettore od una matrice di parametri della verosimiglianza dello
##'   specifico studio.
##' @author Marco Bressan
##' @export
##' @rdname crr-params-norm-low
crr.mean <- function(params, design = c(0, 1)) {
  stopifnot("lunghezza sbagliata per `alpha`" =
              length(params$alpha) == length(design) - (0 %in% design),
            "lunghezza sbagliata per `beta`" =
              length(params$alpha) == length(params$beta))
  if (!0 %in% design) {
    if (USA_MIA_MODELLAZIONE) {
      message("Chiamata mia implementazione!")
      # in questa parte del codice voglio provare ad implementare la mia
      # versione della ncrr senza baseline. si può cambiare settando la
      # variabile globale (a livello di pacchetto) pari a FALSE: in questo caso
      # si userà la parametrizzazione di Guolo !!! SI ASSUME CHE b SIA IN PRIMA
      # POSIZIONE !!!
      mu_ib <- params$alpha[1] + params$beta[1] * params$mu0
      mu <- c(mu_ib, params$alpha[-1] - params$alpha[1] + (params$beta[-1] - params$beta[1]) * mu_ib)
      return(mu)
    }
    # a dispetto del nome, integra anche il calcolo delle medie negli studi
    # non-baseline
    #browser()
    stopifnot(length(params$alpha) == 2) #TODO: se più di 2 studi baseline?!
    return(c(diff(params$alpha), diff(rev(params$alpha))) + c(diff(params$beta), diff(rev(params$beta))) * params$mu0)
  }
  c(0, params$alpha) + c(1, params$beta) * params$mu0
}

##' @rdname crr-params-norm-low
crr.vcov <- function(params, design = c(0, 1)) {
  if (min(design) > 0) {
    if (USA_MIA_MODELLAZIONE) {
      # in questa parte del codice voglio provare ad implementare la mia
      # versione della ncrr senza baseline. si può cambiare settando la
      # variabile globale (a livello di pacchetto) pari a FALSE: in questo caso
      # si userà la parametrizzazione di Guolo !!! SI ASSUME CHE b SIA IN PRIMA
      # POSIZIONE !!!
      betab <- params$beta[-1] - params$beta[1]
      sigmab <- sqrt(params$rho * params$sigma2[-1] * params$sigma2[1])
      if (any(!is.finite(sigmab))) {
        warning("studio baseline fa venire varianze negative!")
        browser()
      }
      vub <- params$beta[1]^2 * params$sigma20 + params$sigma2[1]
      vv <- tcrossprod(c(1, betab)) * vub
      vv[-1, -1] <- vv[-1, -1] + diagoffdiag(1, params$rho, length(design) - 1) * tcrossprod(sigmab)
      return(vv)
    } else {
      #message("Calcolo varcov con studio senza baseline!")
      stopifnot("Covarianza non implementata per design non baseline e più di due studi!" = length(design) <= 2)
    }
    params$beta <- c(params$beta[1] - params$beta[2], params$beta[2] - params$beta[1])
    # sigma_12 = cov(eps1_01, eps1_02) = ???
    return(tcrossprod(params$beta) * params$sigma20 + diag(params$rho * params$sigma2[1] * params$sigma2[2], 2))
  }
  # minore matrice vcov ottenuto togliendo la prima riga e la prima colonna
  #if (length(beta) > 1 && all(c(beta, sigma20, sigma2) != 1)) browser()
  params$sigma2 <- sqrt(params$sigma2)
  #if (any(!is.finite(sigma2))) browser()#stop("sigma2 negativo!")
  vv <- tcrossprod(c(1, params$beta)) * params$sigma20
  vv[-1, -1] <- vv[-1, -1] + tcrossprod(params$sigma2) *
    diagoffdiag(1, params$rho, length(design) - 1)
  vv
}


crr.vcov.achana <- function(sigma20, sigma2, beta, ..., design = c(0, 1)) {
  .Defunct("crr.vcov")
  if (length(params$sigma2) > 1)
    warning("sigma2 ha lunghezza > 1!!")
  params$sigma2 <- rep(params$sigma2[1], length(design) - (0 %in% design))
  params$rho <- 0.5
  crr.vcov(params, design = design)
}


crr.vcov.equivar <- function(sigma20, rho, beta, ..., design = c(0, 1)) {
  .Defunct("crr.vcov")
  if (min(design) > 0) {
    if (USA_MIA_MODELLAZIONE) {
      .NotYetImplemented()
    } else {
      #message("Calcolo varcov con studio senza baseline!")
      stopifnot("Covarianza non implementata per design non baseline e più di due studi!" = length(design) <= 2)
      beta <- c(beta[1] - beta[2], beta[2] - beta[1])
      # sigma_12 = cov(eps1_01, eps1_02) = ???
      return(tcrossprod(beta) * sigma20 + diag(rho * sigma20, 2))
    }
  }
  vv <- tcrossprod(c(1, beta)) * sigma20
  corr <- diag(length(design) - 1) + rho * (1 - diag(length(design) - 1))
  vv[-1, -1] <- vv[-1, -1] + sigma20 * corr
  vv
}


crr.vcov.simple <- function(beta, sigma20, ..., design = c(0, 1)) {
  .Defunct("crr.vcov")
  crr.vcov.equivar(sigma20, .5, beta, design = design)
}


crr.par.idx <- function(np, fixed = NULL, parlen = NULL, lengths = FALSE) {
  plens <- c(alpha = np, beta = np, mu0 = 1, sigma20 = 1, rho = 1, sigma2 = np)
  if (!is.null(parlen)) {
    plens[names(parlen)] <- parlen
  }
  if (length(fixed) > 0)
    plens <- plens[setdiff(names(plens), fixed)]
  ret <- if (lengths) plens else cumsum(plens)
  #cat("DEBUG crr.par.idx: ritorna valori: \n")
  #print(ret)
  return(ret)
}

##' Trasforma i parametri rispettandone il formato (lista o vettore)
##'
##' @title Trasformazione dei parametri
##' @param params Vettore o lista dei parametri
##' @param which Se `params` non ha i nomi, indica di quale parametro si tratta
##'   (assumendo che tutto il vettore appartenga ad uno stesso parametro, ad
##'   esempio `alpha`)
##' @param inverse Effettua la trasformazione inversa (da scala trasformata ad
##'   originale)?
##' @return I parametri trasformati, nello stesso formato di `params` (lista o
##'   vettore).
##' @author Marco Bressan
##' @rdname crr-params
crr.transform.par <- function(params, which = NULL, inverse = FALSE) {
  if (is.null(params))
    return()
  if (is.null(which)) {
    stopifnot("`params` senza nomi" = !is.null(names(params)))
    which <- names(params)
  }
  trans <- .mapply(.ptrans,
                   list(params, names(params)),
                   MoreArgs = list(inverse = inverse))
  names(trans) <- names(params)
  if (!is.list(params))
    trans <- unlist(trans)
  trans
}

.ptrans <- function(x, which, inverse) {
  if (grepl("sigma", which))
    if (inverse) exp(x) else log(x)
  else if (grepl("rho", which))
    if (inverse) (exp(x) - 1) / (exp(x) + 1) else log(1 + x) - log(1 - x)
  else
    x
}
##' @param what Nomi dei parametri da rimuovere
##' @export
##' @rdname crr-params
crr.remove.par <- function(params, what) {
  rem <- c(sigma2 = "sigma2[1-9]+", sigma20 = "sigma20",
           mu0 = "mu0", alpha = "alpha", beta = "beta",
           rho = "rho")
  what <- match.arg(what, names(rem), several.ok = TRUE)
  params[-grep(paste(rem[what], collapse = "|"), names(params))]
}

.parsplit1 <- function(params, pposs, transform = grep("sigma|rho", names(pposs))) {
  plist <- as.list(pposs)
  for (i in seq_along(pposs)) {
    pim1 <- if (i == 1) 1 else pposs[i - 1] + 1
    plist[[i]] <- params[pim1:pposs[i]]
  }
  for (i in transform) {
    plist[[i]] <- .ptrans(plist[[i]], names(plist)[i], TRUE)
  }
  lapply(plist, unname)
}

##' Lo scopo di questa funzione è fornire un tramite tra la rappresentazione
##' dei parametri sottoforma di vettore voluta da `optim` e una più
##' "user-friendly" in cui i parametri sono separati
##'
##' @title Separa i parametri di una NCRR
##' @param params vettore numerico con i parametri in formato "lungo"
##' @param np lunghezza degli alpha - 1. Settando a `NA` si tenta di calcolarlo
##'   sulla base dei nomi del vettore dei parametri.
##' @param transform ritorna i parametri trasformati in scala logaritmica?
##' @param fixed elenco di parametri da escludere. È raccomandato che si
##'   utilizzi l'output di `match.vcov.fixed()`
##' @param parlen Vettore di interi avente i nomi corrispondenti ai parametri e
##'   la lunghezza degli stessi. Se non specificato, legge l'attributo 'parlen'
##'   da `fixed`. ATTENZIONE: per l'utilizzo standard, si consiglia di non
##'   settare esplicitamente questo parametro, ma affidarsi a `fixed`.
##' @author Marco Bressan
##' @export
crr.split.par <- function(params, np = NA, transform = FALSE, fixed = NULL,
                          parlen = attr(fixed, "parlen")) {
  ## if (length(fixed) == 1 && is.na(fixed))
  ##   fixed <- setdiff(PNAMES,
  ##                    stringr::str_extract(names(params),
  ##                                         paste(PNAMES, collapse = "|")))
  if (anyNA(np))
    np <- length(grep("beta", names(params)))
  pposs <- crr.par.idx(np, fixed = fixed, parlen = parlen)
  plist <- NULL
  if (is.matrix(params)) {
    if (max(pposs) != ncol(params)) {
      cat("ERRORE! Stampo gli indici dei parametri che mi aspetto di trovare:")
      print((pposs))
      print(colnames(params))
      if (!is.null(names(params)) && length(ga <- grep("alpha", names(params))) != np)
        warning(sprintf("Set np to %s, but dim(alpha) = %s", np, length(gg)))
      #... funzione per controllo lunghezza parametri
      stop("`params` has wrong dimesion ",
           sprintf("(expected %d, got %d)", max(pposs), ncol(params)))
    }
    extrfn <- \(x, i) x[, i]
  } else {
    if (max(pposs) != length(params)) {
      cat("ERRORE! Stampo gli indici dei parametri che mi aspetto di trovare:")
      print((pposs))
      print(names(params))
      stop("`params` has wrong dimesion ",
           sprintf("(expected %d, got %d)", max(pposs), length(params)))
    }
    extrfn <- \(x, i) x[i]
  }
  for (i in seq_along(pposs)) {
    pim1 <- if (i == 1) 1 else pposs[i - 1] + 1
    plist[[names(pposs)[i]]] <- extrfn(params, pim1:pposs[i])
  }
  if (isTRUE(transform)) {
    for (i in grep("sigma|rho", names(plist), value = TRUE)) {
      plist[[i]] <- .ptrans(plist[[i]], i, TRUE)
    }
  }
  lapply(plist, unname)
}

##' @param ... Può essere usato al posto di `params` per specificare i parametri
##' @export
##' @rdname crr-params
crr.join.par <- function(params, ..., transform = FALSE) {
  if (missing(params)) {
    params <- list(...)[PNAMES]
  }
  params <- lapply(params, unname)
  if (isTRUE(transform)) {
    for (i in c("sigma20", "sigma2", "rho")) {
      params[[i]] <- crr.transform.par(params[[i]], i, FALSE)
    }
  }
  ret <- do.call(if (any(sapply(params, is.matrix))) cbind else c, params)
  drop(ret)
}

fisher.corr <- function() {
  list(
    link = "pearson",
    linkfun = function(mu) 0.5 * log(1 + mu) - log(1 - mu),
    linkinv = function(eta) (exp(2 * eta) - 1) / (exp(2 * eta) + 1)
  )
}

PNAMES <- c("alpha", "beta", "mu0", "sigma20", "rho", "sigma2")
