##' Verosimiglianza per un singolo studio per NCRR __con baseline__.
##' --VERSIONE PRELIMINARE--
##'
##' @title Verosimiglianza NCRR
##' @param params vettore dei parametri contenuti in un unico vettore; in ordine: alfa_0j, beta_0j, mu0, sigma^2_0, rho, sigma^2_{0j}.
##' @param y vettore dei vartheta
##' @param Gamma matrice Gamma
##' @param design vettore che specifica il design (attualmente utile solo la sua dimensione)
##' @return valore della log-verosimiglianza
##' @author Marco Bressan
##' @export
llik1 <- function(params, y, Gamma, design = c(0, 1)) {
  params <- crr.split.par(params, np <- length(design) - 1)
  mu <- rep(NA, np + 1)
  Sigma <- matrix(0, np + 1, np + 1)
  ll <- with(params, {
    mvtnorm::dmvnorm(
      y, mean = mu <<- crr.mean.baseline0(alpha, beta, mu0, design),
      sigma = Sigma <<- crr.vcov(beta, sigma20, sigma2, rho, design) + Gamma,
      log = TRUE)
  })
  #browser()
  cat("PARAMS: ")
  with(params, {
    lapply(list(alpha, beta, mu0, sigma20, rho, sigma2), \(x) deparse1(round(x, 6))) |>
      append(x = list(fmt = "alpha = %s, beta = %s, \n\tmu0 = %s, sigma20 = %s, rho = %s, \n\tsigma2 = %s"), values = _) |>
      do.call(sprintf, args = _) |>
      cat("\nVALUE: ", ll, "\n\n")
    print(list(marginal_mu = mu, marginal_tildeSigma = Sigma))
    cat("\n=====================================\n")
  })
  return(ll)
}

##' Funzione per il calcolo della verosimiglianza per NCRR.
##'
##' @param object oggetto che definisce il design
##' @return una funzione con parametro come `llik1`
##' @export
get.llik.from.design <- function(object, transform = TRUE, echo = 0,
                                 vcov.type = "normal") {
  np <- length(tt <- unique(do.call(c, object$design)))
  fixed.default <- NULL
  if (vcov.type %in% c("simplified", "achana")) {
    fixed.default[["sigma2"]] <- rep(1, np - 1)
  }
  if (vcov.type == "achana") {
    fixed.default[["rho"]] <- 0
  }
  # ritorna la funzione da ottimizzare
  function(params, y = crr.get.theta(object, raw = TRUE),
           Gamma = crr.get.Gamma(object, raw = TRUE),
           fixed = fixed.default) {
    params <- subst.params(crr.split.par(params, np - 1, transform = transform), fixed)

    mu <- do.call(crr.get.mu,
                  append(list(object = object, raw = TRUE),
                         params[c("alpha", "beta", "mu0")]))
    Sigma <- do.call(crr.get.sigma,
                     append(list(object = object, raw = TRUE, type = vcov.type),
                            params[c("beta", "sigma20", "sigma2", "rho")]))
    ll <- sum(mapply(\(t, m, s, g) mvtnorm::dmvnorm(t, mean = m, sigma = s + g,
                                                    log = TRUE),
                     y, mu, Sigma, Gamma))

    if (echo > 0) {
      cat("CURRENT VALUE: ", ll, "\n")
      if (echo > 1) {
        cat("PARAMS: ")
        with(params, {
          lapply(list(alpha, beta, mu0, sigma20, rho, sigma2),
                 \(x) deparse1(round(x, 6))) |>
            append(x = list(fmt = paste("alpha = %s, beta = %s,",
                                        "\tmu0 = %s, sigma20 = %s, rho = %s,",
                                        "\tsigma2 = %s", sep = "\n")), values = _) |>
            do.call(sprintf, args = _) |>
            cat("\n")
        })
        if (echo > 2)
          print(list(marginal_mu = do.call(c, mu),
                     marginal_tildeSigma = blockdiag(Sigma) + blockdiag(Gamma)))
        cat("=====================================\n")
      }
    }
    return(ll)
  }
}
