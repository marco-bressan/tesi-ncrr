### funzioni di prova ########
# studio multi-design con baseline '0'

# alpha = vettore con gli \alpha_bj, idem per beta, mu0 media di \xi
crr.mean.baseline0 <- function(alpha, beta, mu0, design = c(0, 1)) {
  stopifnot("wrong size for alpha" = length(alpha) == length(design) - 1,
            "wrong size for beta" = length(alpha) == length(beta))
  c(0, alpha) + c(1, beta) * mu0
}

# sigma2 vettore varianze di dim compatibile con beta
crr.vcov.baseline0 <- function(beta, sigma20, sigma2, rho, design = c(0, 1)) {
  # minore matrice vcov ottenuto togliendo la prima riga e la prima colonna
  vv1sigma2 <- rho * outer(sigma2, sigma2)
  diag(vv1sigma2) <- 0
  vv1 <- outer(beta, beta) * (sigma20 + vv1sigma2)
  diag(vv1) <- diag(vv1) + sigma2
  # costruzione output
  vv <- matrix(NA, length(design), length(design))
  vv[1, ] <- c(1, beta) * sigma20
  vv[-1, 1] <- vv[1, -1]
  vv[-1, -1] <- vv1

  vv
}

crr.vcov.within <- function(r, n) {
  diag(1 / r + 1 / (n + r))
}

crr.split.par <- function(params, np) {
  params <- unname(params)
  alpha <- params[1:np]
  beta <- params[(1 + np):(2 * np)]
  mu0 <- params[2 * np + 1]
  sigma20 <- params[2 * np + 2]
  rho <- params[2 * np + 3]
  sigma2 <- params[(2 * np + 4):length(params)]
  #print(as.list(environment()))
  stopifnot("`params` has wrong dimesion!" = length(sigma2) == np)
  rm(np, params)
  return(as.list(environment()))
}
