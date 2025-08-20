random.design <- function() {

}

#' Funzione per simulare da un design. Può essere usato in combinazione con
#' `random.design()` per svolgere studi di simulazione approfonditi.
#'
#' @param object Un oggetto di classe `ncrr.design`
#' @param nsim numero di simulazioni
#' @param seed Seme casuale per riproducibilità
#' @param vcov.type Tipo di matrice di varianza-covarianza da usare
#' @param ... Parametri previsti per la specifica tipologia di matrice
#' (un sottoinsieme di alpha, beta, mu0, sigma20, rho, sigma2)
#'
#' @returns Una lista di oggetti `ncrr.design` con i dati simulati.
#' @export
#'
simulate.ncrr.design <- function(object, nsim = 1, seed = rpois(1, 1e5),
                                 vcov.type = attr(object, "vcov.type"),
                                 params, ...) {
  if (is.null(vcov.type))
    vcov.type <- "normal"
  vcov.type <- match.vcov.type(vcov.type)
  set.seed(seed)
  dd <- object$design
  ds <- object$x
  if (missing(params))
    params <- list(...)
  params <- set.vcov.params(params, vcov.type = vcov.type)
  mu <- crr.get.mu(object, params, raw = TRUE)
  #browser()
  Sigma <- crr.get.sigma(object, params, raw = TRUE)
  theta <- mapply(\(m, S) t(rmvnorm(nsim, m, S)), mu, Sigma, SIMPLIFY = FALSE)
  theta <- do.call(rbind, theta)
  pik <- exp(theta) / (1 + exp(theta))
  rik <- matrix(rbinom(length(pik), rep(ds$nik, ncol(pik)), c(pik)),
                nrow(pik), ncol(pik))
  structure(
    lapply(seq_len(ncol(pik)), \(i) {
      ds2 <- ds
      ds2$tik <- theta[, i]
      ds2$rik <- rik[, i]
      ncrr.design(ds2, vcov.type)
    }),
    seed = seed,
    class = "ncrr.simul")
}

#' @export
print.ncrr.simul <- function(x, ...) {
  seed <- attr(x, "seed")
  cat("Simulazione NCRR", if (!is.null(seed)) sprintf("(seed = %i)", seed))
  cat("avente ", length(x), " replicazion", ifelse(length(x) == 1, "e", "i"),
      ".\n", sep = "")
  print(x[[1]])
  cat("\n")
}
