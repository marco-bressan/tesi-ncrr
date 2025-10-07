##' Calcolo della statistica di Skovgaard
##'
##'
##' @title Calcolo della statistica di Skovgaard
##' @param data Un oggetto di classe `ncrr.design`
##' @param thetainit Stima iniziale del vettore dei parametri
##' @param floglik Funzione di verosimiglianza, idealmente il prodotto di una
##'   chiamata a `get.llik.from.design()`
##' @param fscore Funzione score, se disponibile (altrimenti viene usata
##'   l'approssimazione numerica del gradiente contenuta in `pracma`)
##' @param fpsi Funzione scalare che estrapola il valore del parametro
##'   d'interesse
##' @param psival Valore del parametro d'interesse sotto l'ipotesi nulla
##' @param datagen Funzione che genera il dataset a partire dalla stima
##'   corrente dei parametri
##' @param R Numero di iterazioni per il calcolo della statistica di Skovgaard
##'   tramite Montecarlo
##' @param seed Seed per il metodo Montecarlo
##' @param ronly Se TRUE, non calcola la statistica di Skovgaard
##' @param parallel Calcolo parallelo per la statistica di Skovgaard
##' @param nclus Numero di cluster per il calcolo parallelo
##' @param trace Intero >= 0
##' @param precalc Lista di valori precalcolati (SPERIMENTALE!!)
##' @return *Da completare...*
##' @author Marco Bressan
##' @export
crr.rstar <- function(data, thetainit, floglik, fscore = NULL, fpsi, psival,
                      datagen, R = 1000, seed = NULL, ronly = FALSE,
                      parallel = FALSE, nclus = NA, trace = 1, precalc = NULL) {
  if (!is.list(data)) {
    warning("data should be provided as a list\n")
    data <- as.list(data)
  }
  if (!is.numeric(thetainit))
    stop("a starting point for the parameter theta is required \n")
  if (!is.numeric(psival))
    stop("a value for the parameter of interest is required \n")
  f0 <- floglik(thetainit, data)
  if (!is.numeric(f0))
    stop("problems in the loglikelihood function \n")
  if (!is.null(fscore)) {
    g0 <- fscore(thetainit, data)
    if (length(g0) != length(thetainit))
      stop("size of starting point different from the size of score function \n")
  }
  if (!ronly && !is.numeric(floglik(thetainit, datagen(data, thetainit)))) {
    stop("problems in the function to simulate data \n")
  }
  p <- length(thetainit)
  if (trace > 0)
    message("get mle ....", "\t")
  min.floglik <- \(theta, data) -floglik(theta, data)
  min.fscore <- if (!is.null(fscore)) \(theta, data) -fscore(theta, data)
  el.hat <- -Inf
  el.til <- Inf
  count <- 0
  while (count > 100 || el.hat < el.til) {# parziale workaround per minimi locali
    obj.hat <- nlminb(thetainit, min.floglik, min.fscore, data = data)
    theta.hat <- obj.hat$par
    el.hat <- floglik(theta.hat, data)
    objHyp <- Rsolnp::solnp(theta.hat, fun = min.floglik,
                            eqfun = \(theta, data) fpsi(theta), eqB = psival,
                            control = list(trace = 0), data = data)
    theta.til <- objHyp$par
    el.til <- floglik(theta.til, data)
    thetainit <- theta.til
    count <- count + 1
  }
  if (trace > 0)
    message("mle attained in ", count, " iteration(s)")

  psi.hat <- fpsi(theta.hat)
  r <- sqrt(2 * (el.hat - el.til)) * sign(psi.hat - psival)
  out <- list(r = r, theta.hat = theta.hat, psi.hat = psi.hat,
              theta.hyp = theta.til, psi.hyp = psival)
  class(out) <- "crr.rstar"
  if (ronly) return(out)

  if (abs(r) < 0.1)
    warning("Value under testing close to the MLE - ",
            "there might be a singularity in r*\n")

  # calcolo matrici di informazione osservata e varianze
  j.hat <- {
    if (is.null(fscore))
      -pracma::hessian(floglik, theta.hat, data = data)
    else
      -pracma::jacobian(fscore, theta.hat, data = data)
  }
  score.hat.data <- {
    if (is.null(fscore))
      pracma::grad(floglik, theta.hat, data = data)
    else
      fscore(theta.hat, data = data)
  }
  var.theta.hat <- solve(j.hat)
  se.theta.hat <- sqrt(diag(var.theta.hat))
  j.til <- {
    if (is.null(fscore))
      -pracma::hessian(floglik, theta.til, data = data)
    else
      -pracma::jacobian(fscore, theta.til, data = data)
  }
  score.til.data <- {
    if (is.null(fscore))
      pracma::grad(floglik, theta.til, data = data)
    else
      fscore(theta.til, data)
  }
  dpsi.dtheta <- pracma::grad(fpsi, theta.hat)
  var.psi.hat <- dpsi.dtheta %*% var.theta.hat %*%
    dpsi.dtheta
  se.psi.hat <- sqrt(var.psi.hat)
  C.hat <- pracma::grad(fpsi, theta.hat)
  C.til <- pracma::grad(fpsi, theta.til)
  k <- which(C.hat != 0)[1]
  obj.info <- likelihoodAsy:::.newinfo(p, k, C.hat, C.til, j.hat,
                                       j.til, score.hat.data, score.til.data, theta.hat,
                                       theta.til, fpsi)
  j.hat.new <- obj.info$j.hat.new
  j.til.new <- obj.info$j.til.new
  # calcolo rstar tramite montecarlo
  input <- .mapply(\(...) list(...), append(list(1:R), precalc), NULL)
  Uhh <- .parallel(
    parallel = parallel, nclus = nclus, trace = trace, seed = seed,
    initfun = \() devtools::load_all(),
    pb = txtProgressBar(style = 3),
    trace.init.msg = "Starting Monte-carlo simulation",
    x = input,
    fun = \(x) {
      dataSim <- x$data %||% datagen(data, theta.hat)
      l1 <- floglik(theta.hat, dataSim)
      l0 <- floglik(theta.til, dataSim)
      score.hat <- {
        if (is.null(fscore))
          pracma::grad(f = floglik, x0 = theta.hat, data = dataSim)
        else
          fscore(theta.hat, dataSim)
      }
      score.til <- {
        if (is.null(fscore))
          pracma::grad(f = floglik, x0 = theta.til, data = dataSim)
        else
          fscore(theta.til, dataSim)
      }
      obj.score <- likelihoodAsy:::.newscores(p, k, C.hat, C.til,
                                              score.hat, score.til)
      c(obj.score$score.new.hat, obj.score$score.new.til, l1 - l0)
    }
  )
  meanAll <- Reduce("+", Uhh, init = rep(0, 2 * p + 1)) / R
  prodAll <- Reduce(\(u1, u2) u1 + tcrossprod(u2), Uhh,
                    init = matrix(0, 2 * p + 1, 2 * p + 1)) / R
  covAll <- prodAll * R/(R - 1) - tcrossprod(meanAll) * R/(R - 1)
  S <- covAll[1:p, (p + 1):(2 * p)]
  i.hat <- covAll[1:p, 1:p]
  i.hatInv <- qr.solve(i.hat, tol = 10^-20)
  q <- covAll[2 * p + 1, 1:p]
  SS2 <- t(S) %*% i.hatInv %*% j.hat.new
  SS1 <- q %*% i.hatInv %*% j.hat.new
  indpsi <- k
  numU <- SS1[indpsi] -
    SS2[-indpsi, indpsi] %*% solve(t(SS2[-indpsi, -indpsi])) %*% SS1[-indpsi]
  j.hatInv.new <- solve(j.hat.new)
  jProf.new <- 1/j.hatInv.new[indpsi, indpsi]
  u <- numU/sqrt(jProf.new)
  CPsi <- det(as.matrix(SS2[-indpsi, -indpsi])) /
    sqrt(det(as.matrix(j.til.new[-indpsi, -indpsi])) *
           det(as.matrix(j.hat.new[-indpsi, -indpsi])))
  NP <- (1/r) * log(CPsi)
  INF <- (1/r) * log(u/r)
  rs <- r + NP + INF
  out$NP <- drop(NP)
  out$INF <- drop(INF)
  out$rs <- drop(rs)
  out$info.hat <- j.hat
  out$se.theta.hat <- se.theta.hat
  out$se.psi.hat <- drop(se.psi.hat)
  out$seed <- attr(Uhh, "seed")
  out$R <- R

  out
}
