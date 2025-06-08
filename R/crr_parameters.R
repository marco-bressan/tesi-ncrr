### funzioni di prova ########
# studio multi-design con baseline '0'

# alpha = vettore con gli \alpha_bj, idem per beta, mu0 media di \xi
crr.mean.baseline0 <- function(alpha, beta, mu0, design = c(0, 1)) {
  stopifnot("wrong size for alpha" = length(alpha) == length(design) - 1,
            "wrong size for beta" = length(alpha) == length(beta))
  c(0, alpha) + c(1, beta) * mu0
}

# sigma2 vettore varianze di dim compatibile con beta
crr.vcov <- function(beta, sigma20, sigma2, rho, design = c(0, 1)) {
  if (min(design) > 0) {
    stop("Covarianza non implementata per design non baseline!")
  }
  # minore matrice vcov ottenuto togliendo la prima riga e la prima colonna
  #if (length(beta) > 1 && all(c(beta, sigma20, sigma2) != 1)) browser()
  sigma2 <- sqrt(sigma2)
  if (any(!is.finite(sigma2))) stop("sigma2 negativo!")
  vv <- tcrossprod(c(1, beta)) * sigma20
  vv[-1, -1] <- vv[-1, -1] + tcrossprod(sigma2) *
    (rho * (1 - diag(length(design) - 1)) + diag(length(design) - 1))

  vv
}

crr.vcov.simple <- function(beta, sigma20, rho, ..., design = c(0, 1)) {
  ...not.used(...)
  if (min(design) > 0) {
    stop("Covarianza non implementata per design non baseline!")
  }
  vv <- tcrossprod(c(1, beta)) * sigma20
  vv[-1, -1] <- vv[-1, -1] + sigma20 *
    (rho * (1 - diag(length(design) - 1)) + diag(length(design) - 1))

  vv
}

crr.vcov.achana <- function(sigma20, ..., design = c(0, 1)) {
  ...not.used(...)
  stopifnot(length(sigma20) == 1)
  sigma20 * (diag(length(design)) + .5 * (1 - diag(length(design))))
}

crr.vcov.within <- function(r, n) {
  diag(1 / r + 1 / (n + r))
}

crr.split.par <- function(params, np, transform = FALSE) {
  if (is.matrix(params)) {
    alpha <- params[, 1:np]
    beta <- params[, (1 + np):(2 * np)]
    mu0 <- params[, 2 * np + 1]
    sigma20 <- params[, 2 * np + 2]
    rho <- params[, 2 * np + 3]
    sigma2 <- params[, (2 * np + 4):ncol(params)]
    stopifnot("`params` has wrong dimesion!" = ncol(sigma2) == np)
  } else {
    params <- unname(params)
    alpha <- params[1:np]
    beta <- params[(1 + np):(2 * np)]
    mu0 <- params[2 * np + 1]
    sigma20 <- params[2 * np + 2]
    rho <- params[2 * np + 3]
    sigma2 <- params[(2 * np + 4):length(params)]
    stopifnot("`params` has wrong dimesion!" = length(sigma2) == np)
  }

  #print(as.list(environment()))

  if (isTRUE(transform)) {
    sigma20 <- exp(sigma20)
    sigma2 <- exp(sigma2)
    rho <- fisher.corr()$linkinv(rho)
  }
  rm(np, params, transform)
  return(as.list(environment()))
}

crr.join.par <- function(params, ..., transform = FALSE) {
  pnames <- c("alpha", "beta", "mu0", "sigma20", "rho", "sigma2")
  if (!missing(params)) {
    stopifnot("alcuni nomi mancanti in `params`" = pnames %in% names(params))
  } else {
    if (any(!pnames %in% ...names()))
      stop("Non tutti i parametri sono stati passati correttamente!")
    params <- list(...)[pnames]
  }
  if (isTRUE(transform)) {
    params[["sigma20"]] <- log(params[["sigma20"]])
    params[["sigma2"]] <- log(params[["sigma2"]])
    params[["rho"]] <- fisher.corr()$linkfun(params[["sigma20"]])
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
