#' Questa funzione calcola la statistica radice con segno ($r_p$), basata
#' sul log-rapporto di verosimiglianza, utile nell'inferenza parametrica
#' condizionale e nell'analisi della significatività di parametri di interesse.
#'
#' @title Calcolo della statistica log-RV
#'
#' @param dati.gen Oggetto di classe ncrr.design
#' @param psi0 Valore del parametro di interesse $\psi$ sotto l'ipotesi nulla.
#' @param init Vettore di valori iniziali per l'ottimizzazione dei parametri.
#' @param param Indice del parametro di interesse $\psi$. Di default si
#'   considera il parametro "beta5".
#' @param theta.hat (Opzionale) Punto di massima della verosimiglianza
#'   precalcolato. Se NULL viene stimata internamente.
#' @param J (Opzionale) Matrice hessiana della funzione di log-verosimiglianza
#'   precalcolata. Se NULL viene stimata internamente.
#' @param ... Argomenti aggiuntivi passati alla funzione di verosimiglianza.
#' @param exact Se TRUE, usa ottimizzazione vincolata esatta.
#' @param par.only Se TRUE, restituisce solo i parametri ottimizzati sotto
#'   vincolo
#' @param export Nomi degli attributi da esportare: di default, `theta.hat`
#'   (stima libera dei parametri), `l.hat` (verosimiglianza nel punto di
#'   massimo) e `J` (matrice hessiana)
#' @return La statistica $r_p$ come oggetto numerico con gli attributi
#'   specificati in `export`. Se `par.only = TRUE`, restituisce il vettore dei
#'   parametri ottimizzati sotto vincolo.
#' @export
rp.stat <- function(dati.gen, psi0, init, param = match("beta5", names(init)),
                    theta.hat = NULL, J = NULL, ...,
                    exact = NA, par.only = FALSE, export = c("theta.hat", "l.hat", "J")) {
  # psi par d'interesse, lam di disturbo
  if (is.null(theta.hat)) {
    opt.theta <- optim(init, \(x) -llik.fun(x, dati.gen), method = "BFGS", hessian = TRUE)
    l.hat <- -opt.theta$value
    theta.hat <- opt.theta$par
    J <- opt.theta$hessian
  } else {
    l.hat <- llik.fun(theta.hat, dati.gen)
    if (is.null(J))
      J <- -optimHess(theta.hat, llik.fun)
  }
  # opt.theta <- optim(init, \(x) -llik.fun(x, dati.gen), method = "BFGS", hessian = TRUE)
  # if (!isTRUE(all.equal(opt.theta$par, theta.hat)) ||
  #     !isTRUE(all.equal(opt.theta$hessian, J)))
  #   browser()
  theta.psi <- theta.hat
  if (!isTRUE(exact)) {
    theta.psi[param] <- psi0
    # calcolo approssimato del parametro di disturbo
    theta.psi[-param] <- theta.hat[-param] +
      c(solve(J[-param, -param]) %*% J[-param, param] %*% (theta.hat[param] - psi0))
    ## confronto
    # theta.psi2 <- Rsolnp::solnp(theta.hat, \(x) -llik.fun(x, dati.gen),
    #                            eqfun = \(t) t[param], eqB = psi0,
    #                            control = list(trace = 0))$pars
    # print(paste("Psi0 =", psi0))
    # print(all.equal(theta.psi, theta.psi2))
    # print(((theta.psi - theta.psi2)/theta.psi))
  }
  if (!isFALSE(exact)) {
    theta.psi <- Rsolnp::solnp(theta.psi, \(x) -llik.fun(x, dati.gen),
                               eqfun = \(t) t[param], eqB = psi0,
                               control = list(trace = 0))$pars
    ## Equivalente ma più lenta:
    # lam0 <- nlminb(theta.hat[-param], \(x) {
    #   z <- theta.hat
    #   z[param] <- psi0
    #   z[-param] <- x
    #   -llik.fun(z, dati.gen)
    # })$par
    # theta.psi <- theta.hat
    # theta.psi[param] <- psi0
    # theta.psi[-param] <- lam0
  }
  if (par.only) return(theta.psi)
  lp0 <- llik.fun(theta.psi, dati.gen)
  rp <- unname(sign(theta.hat[param] - psi0) * sqrt(2) * sqrt(l.hat - lp0))
  for (a in export) attr(rp, a) <- get(a, environment())
  rp
}
