crr.rstar <- function(data, thetainit, floglik, fscore = NULL, fpsi, psival,
                      datagen, R = 1000, seed = NULL, ronly = FALSE,
                      psidesc = NULL, constr.opt = "solnp",
                      parallel = FALSE, nclus = NA, trace = 1) {
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
  if (!ronly) {
    data0 <- datagen(thetainit, data)
    f0 <- floglik(thetainit, data0)
    if (!is.numeric(f0))
      stop("problems in the function to simulate data \n")
  }
  if ((constr.opt != "solnp") & (constr.opt != "alabama"))
    stop("constrained optimizer must be either 'solnp' or 'alabama'")
  p <- length(thetainit)
  if (trace > 0)
    cat("get mle ....", "\t")
  min.floglik <- \(theta, data) -floglik(theta, data)
  min.fscore <- if (!is.null(fscore)) \(theta, data) -fscore(theta, data)
  obj.hat <- nlminb(thetainit, min.floglik, min.fscore, data = data)
  theta.hat <- obj.hat$par
  el.hat <- floglik(theta.hat, data)
  if (!ronly) {
    j.hat <- if (is.null(fscore))
               -pracma::hessian(floglik, theta.hat, data = data)
    else -pracma::jacobian(fscore, theta.hat, data = data)
    score.hat.data <- if (is.null(fscore))
                        pracma::grad(floglik, theta.hat, data = data)
    else fscore(theta.hat, data = data)
    var.theta.hat <- solve(j.hat)
    se.theta.hat <- sqrt(diag(var.theta.hat))
  }
  if (trace > 0)
    cat("get mle under the null....", "\n")
  psifcn.mod <- if (constr.opt == "solnp")
                  function(theta, data) fpsi(theta)
  else function(theta, data) fpsi(theta) - psival
  objHyp <- if (constr.opt == "solnp")
              Rsolnp::solnp(theta.hat, fun = min.floglik, eqfun = psifcn.mod,
                            eqB = psival, control = list(trace = 0), data = data)
  else alabama::constrOptim.nl(theta.hat, fn = min.floglik, heq = psifcn.mod,
                               gr = min.fscore, control.outer = list(trace = FALSE),
                               data = data)
  theta.til <- objHyp$par
  el.til <- floglik(theta.til, data)
  psi.hat <- fpsi(theta.hat)
  if (!ronly) {
    j.til <- if (is.null(fscore))
               -pracma::hessian(floglik, theta.til, data = data)
    else -pracma::jacobian(fscore, theta.til, data = data)
    score.til.data <- if (is.null(fscore))
                        pracma::grad(floglik, theta.til, data = data)
    else fscore(theta.til, data)
    dpsi.dtheta <- pracma::grad(fpsi, theta.hat)
    var.psi.hat <- dpsi.dtheta %*% var.theta.hat %*%
      dpsi.dtheta
    se.psi.hat <- sqrt(var.psi.hat)
  }
  r <- sqrt(2 * (el.hat - el.til)) * sign(psi.hat - psival)
  if (!ronly) {
    C.hat <- pracma::grad(fpsi, theta.hat)
    C.til <- pracma::grad(fpsi, theta.til)
    k <- which(C.hat != 0)[1]
    obj.info <- likelihoodAsy:::.newinfo(p, k, C.hat, C.til, j.hat,
                                         j.til, score.hat.data, score.til.data, theta.hat,
                                         theta.til, fpsi)
    j.hat.new <- obj.info$j.hat.new
    j.til.new <- obj.info$j.til.new
  }
  if (ronly)
    out <- list(r = r, theta.hat = theta.hat, psi.hat = psi.hat,
                theta.hyp = theta.til, psi.hyp = psival)
  if (!ronly) {
    para <- .setup.parallel(parallel, nclus, trace, seed, R,
                            pb = txtProgressBar(style = 3),
                            trace.init.msg = "Starting Monte-carlo simulation")
    Uhh <- snowFT::performParallel(
      count = para$nclus,
      x = 1:R,
      fun = \(i) {
        dataSim <- datagen(theta.hat, data = data)
        l1 <- floglik(theta.hat, dataSim)
        l0 <- floglik(theta.til, dataSim)
        score.hat <- if (is.null(fscore))
                       pracma::grad(f = floglik, x0 = theta.hat,
                                    data = dataSim)
        else fscore(theta.hat, dataSim)
        score.til <- if (is.null(fscore))
                       pracma::grad(f = floglik, x0 = theta.til, data = dataSim)
        else fscore(theta.til, dataSim)
        obj.score <- likelihoodAsy:::.newscores(p, k, C.hat, C.til,
                                                score.hat, score.til)
        c(obj.score$score.new.hat, obj.score$score.new.til, l1 - l0)
      },
      printfun = para$printfun,
      printrepl = para$printrepl,
      ft_verbose = trace > 1,
      seed = para$seed
    )
    close.parallel(para)
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
    out <- list(r = r, NP = drop(NP), INF = drop(INF),
                rs = drop(rs), theta.hat = theta.hat, info.hat = j.hat,
                se.theta.hat = se.theta.hat, psi.hat = psi.hat,
                se.psi.hat = drop(se.psi.hat), theta.hyp = theta.til,
                psi.hyp = psival, seed = para$seed)
  }
  out$psidesc <- psidesc
  out$R <- R
  if ((!ronly) & (abs(r) < 0.1)) {
    cat("Value under testing close to the MLE - there might be a singularity in r*\n")
    warning("Value under testing close to the MLE - there might be a singularity in r*\n")
  }
  return(structure(out, class = "rstar"))
}
