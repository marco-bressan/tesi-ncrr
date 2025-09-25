#' Calcolo della statistica radice con segno del log-rapporto di verosimiglianza (\( r_p \))
#'
#' Questa funzione calcola la statistica radice con segno (\( r_p \)), basata sul log-rapporto di verosimiglianza, utile nell'inferenza parametrica condizionale e nell'analisi della significatività di parametri di interesse.
#'
#' ## Dettagli teorici
#' \[
#' r_p = \mathrm{sign}(\hat{\psi} - \psi_0) \sqrt{2(\ell(\hat{\theta}) - \ell(\tilde{\theta}))}
#' \]
#' dove \(\ell(\hat{\theta})\) è la log-verosimiglianza nell'ipotesi libera, \(\ell(\tilde{\theta})\) è la log-verosimiglianza vincolata sotto \(\psi = \psi_0\), \(\hat{\psi}\) è la stima massima di \(\psi\) e \(\psi_0\) il valore fissato.
#'
#' ## Argomenti
#' @param dati.gen Lista o vettore di dati osservati da analizzare.
#' @param psi0 Valore del parametro di interesse \(\psi\) sotto ipotesi nulla.
#' @param init Vettore di valori iniziali per l'ottimizzazione dei parametri.
#' @param param Indice (o nome) del parametro di interesse \(\psi\). Default: `match("beta5", names(init))`.
#' @param theta.hat (Opzionale) Stima libera massima della verosimiglianza. Se NULL viene stimata internamente.
#' @param J (Opzionale) Matrice hessiana della funzione di log-verosimiglianza. Se NULL viene stimata internamente.
#' @param ... Argomenti aggiuntivi passati alla funzione di verosimiglianza.
#' @param exact Se TRUE, usa ottimizzazione vincolata esatta; se FALSE, usa approssimazione di Taylor per i disturbi.
#' @param par.only Se TRUE, restituisce solo i parametri ottimizzati sotto vincolo. Default: FALSE.
#'
#' ## Valore
#' Restituisce la statistica \( r_p \) come oggetto numerico, con attributi `theta.hat` (stima libera dei parametri) e `J` (matrice hessiana). Se `par.only = TRUE`, restituisce il vettore dei parametri ottimizzati sotto vincolo.
#'
#' ## Esempio
#' ```
#' # Suppose dati.gen is a vector of outcomes, init is a named vector of initial parameter guesses
#' rp <- rp.stat(dati.gen, psi0 = 0.2, init)
#' ```
#'
#' @export
rp.stat <- function(dati.gen, psi0, init, param = match("beta5", names(init)),
                    theta.hat = NULL, J = NULL, ..., exact = NA, par.only = FALSE) {
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
  structure(rp, theta.hat = theta.hat, J = J)
}
