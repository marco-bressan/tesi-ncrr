### funzioni di prova ########
# studio multi-design con baseline '0'

# alpha = vettore con gli \alpha_bj, idem per beta, mu0 media di \xi
crr.mean.baseline0 <- function(alpha, beta, mu0, design = c(0, 1)) {
  stopifnot("wrong size for alpha" =
              length(alpha) == length(design) - (0 %in% design),
            "wrong size for beta" = length(alpha) == length(beta))
  if (!0 %in% design) {
    # a dispetto del nome, integra anche il calcolo delle medie negli studi non-baseline
    #browser()
    stopifnot(length(alpha) == 2) #TODO: se più di 2 studi baseline?!
    return(c(diff(alpha), diff(rev(alpha))) + c(diff(beta), diff(rev(beta))) * mu0)
  }
  c(0, alpha) + c(1, beta) * mu0
}

# sigma2 vettore varianze di dim compatibile con beta
crr.vcov <- function(beta, sigma20, sigma2, rho, design = c(0, 1)) {
  if (min(design) > 0) {
    warning("Covarianza non implementata per design non baseline!")
    stopifnot("Covarianza non implementata per design non baseline e più di due studi!" = length(design > 2))
    beta <- c(beta[1] - beta[2], beta[2] - beta[1])
    # sigma_12 = cov(eps1_01, eps1_02) = ???
    return(tcrossprod(beta) * sigma20 + diag(rho * sigma2[1] * sigma2[2], 2))
  }
  # minore matrice vcov ottenuto togliendo la prima riga e la prima colonna
  #if (length(beta) > 1 && all(c(beta, sigma20, sigma2) != 1)) browser()
  sigma2 <- sqrt(sigma2)
  #if (any(!is.finite(sigma2))) browser()#stop("sigma2 negativo!")
  vv <- tcrossprod(c(1, beta)) * sigma20
  vv[-1, -1] <- vv[-1, -1] + tcrossprod(sigma2) *
    (rho * (1 - diag(length(design) - 1)) + diag(length(design) - 1))
  vv
}

crr.vcov.simple <- function(beta, sigma20, rho, ..., design = c(0, 1)) {
  #...not.used(...)
  if (min(design) > 0) {
    stopifnot("Covarianza non implementata per design non baseline e più di due studi!" = length(design > 2))
    beta <- c(beta[1] - beta[2], beta[2] - beta[1])
    return(tcrossprod(beta) * sigma20 + diag(sigma20 / 2, 2))
  }
  vv <- tcrossprod(c(1, beta)) * sigma20
  vv[-1, -1] <- vv[-1, -1] + sigma20 *
    (rho * (1 - diag(length(design) - 1)) + diag(length(design) - 1))

  vv
}

crr.vcov.achana <- function(sigma20, ..., design = c(0, 1)) {
  #...not.used(...)
  stopifnot(length(sigma20) == 1)
  #if (runif(1) < .05) browser()
  sigma20 * (diag(length(design)) + .5 * (1 - diag(length(design))))
}

crr.vcov.within <- function(r, n) {
  diag(1 / r + 1 / (n + r))
}

crr.par.idx <- function(np, fixed = NULL, lengths = FALSE) {
  plens <- c(alpha = np, beta = np, mu0 = 1, sigma20 = 1, rho = 1, sigma2 = np)
  if (length(fixed) > 0)
    plens <- plens[setdiff(names(plens), fixed)]
  return(if (lengths) plens else cumsum(plens))
}

crr.transform.par <- function(value, which, inverse = FALSE) {
  if (is.null(value))
    return()
  tfun <- {
    if (grepl("sigma", which)) {
      if (inverse) exp else log
    } else if (grepl("rho", which)) {
      fisher.corr()[[if (inverse) "linkinv" else "linkfun"]]
    } else {
      \(x) x
    }
  }
  tfun(value)
}

crr.remove.par <- function(params, what) {
  rem <- c(sigma2 = "sigma2[1-9]+", sigma20 = "sigma20",
           mu0 = "mu0", alpha = "alpha", beta = "beta",
           rho = "rho")
  what <- match.arg(what, names(rem), several.ok = TRUE)
  params[-grep(paste(rem[what], collapse = "|"), names(params))]
}

##' lo scopo di questa funzione è fornire un tramite tra la rappresentazione
##' dei parametri sottoforma di vettore voluta da `optim` e una più
##' "user-friendly" in cui i parametri sono separati
##'
##' @title separa i parametri di una ncrr
##' @param params vettore numerico con i parametri in formato "lungo"
##' @param np lunghezza degli alpha - 1
##' @param transform ritorna i parametri trasformati in scala logaritmica?
##' @param fixed elenco di parametri da escludere
##' @return una lista con i parametri separati
##' @author Marco Bressan
crr.split.par <- function(params, np, transform = FALSE, fixed = NULL) {
  ## if (length(fixed) == 1 && is.na(fixed))
  ##   fixed <- setdiff(PNAMES,
  ##                    stringr::str_extract(names(params),
  ##                                         paste(PNAMES, collapse = "|")))
  pposs <- crr.par.idx(np, fixed = fixed)
  plist <- NULL
  if (is.matrix(params)) {
    if (max(pposs) != ncol(params)) {
      print((pposs))
      print(colnames(params))
      stop("`params` has wrong dimesion!")
    }
    extrfn <- \(x, i) x[, i]
  } else {
    if (max(pposs) != length(params)) {
      print((pposs))
      print(names(params))
      stop("`params` has wrong dimesion!")
    }
    extrfn <- \(x, i) x[i]
  }
  for (i in seq_along(pposs)) {
    pim1 <- if (i == 1) 1 else pposs[i - 1] + 1
    plist[[names(pposs)[i]]] <- extrfn(params, pim1:pposs[i])
  }
  if (isTRUE(transform)) {
    for (i in c("sigma20", "sigma2", "rho")) {
      plist[[i]] <- crr.transform.par(plist[[i]], i, TRUE)
    }
  }
  lapply(plist, unname)
}

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
    linkfun = function(mu) log(1 + mu) - log(1 - mu),
    linkinv = function(eta) (exp(eta) - 1) / (exp(eta) + 1)
  )
}

PNAMES <- c("alpha", "beta", "mu0", "sigma20", "rho", "sigma2")
