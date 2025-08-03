#' ---
#' title: studio di simulazione
#' subtitle: dati di Achana
#' format:
#'   pdf:
#'     fig-width: 12
#'     fig-height: 10
#' ---
rm(list = ls())
#setwd("/home/marco/UNIPD Scienze Statistiche/tesi-ncrr/data-raw")
devtools::load_all()

CONFLVL <- .95
NSIM <- 250
VCOVTYPE <- "achana"
DIR <- "../.." # per il markdown
#DIR <- ".." # per l'esecuzione nel pacchetto

#' # Simulazione basata sul problema di achana
#| warning: false
simu.pars <- list(alpha = c(0.53118984013899, 1.0431973777787, 0.00434231242384523,
                            2.36407165618289, 2.66293182986318, 2.7339581049579),
                  beta = c(0.948918313002282,
                           1.02425107553256, 1.06604309779969, 0.240340161234799, 0.179383944934584,
                           0.179378353526596),
                  mu0 = 0.81098898311333,
                  sigma20 = 2.63212049308308,
                  sigma2 = 5.69469982077026) # stime MV dai dati originali
simu.des <- do.call(simulate,
                    append(list(ncrr.design(smoke.alarm), vcov.type = VCOVTYPE,
                                nsim = NSIM, seed = 212),
                           simu.pars))
llik.fun <- function(theta, data) {
  y = crr.get.theta(data, raw = TRUE)
  Gamma = crr.get.Gamma(data, raw = TRUE)
  L <- llik(theta, y = y, Gamma = Gamma)
  # if (!is.finite(L)) {
  #   environment(llik)$echo <- 3
  #   browser()
  #   llik(theta, y = y, Gamma = Gamma)
  # }
  return(L)
}
gendat.fun <- function(theta, data) {
  theta <- crr.split.par(theta, length(data$treatments) - 1,
                         fixed = match.vcov.fixed(attr(data, "vcov.type")))
  suppressMessages(simulate(data, params = theta, seed = 342)[[1]]) # estraggo solo il design, non tutta la lista
}
psi.fun <- function(theta) {
  theta[["beta5"]]
}

#| eval: false
# simulazione
lik.vals <- lik.vals2 <- numeric(NSIM)
init <- getInitial(simu.des[[1]], vcov.type = VCOVTYPE)
par.J <- par.J2 <- array(NA, c(length(init), length(init), NSIM))
par.stime <- par.stime2 <- par.sd <- par.sd2 <- matrix(NA, length(init), NSIM)
for (k in seq_len(NSIM)) {
  cat(sprintf("%.2f%%\r", k / NSIM * 100))
  llik <- get.llik.from.design(simu.des[[k]], vcov.type = VCOVTYPE, echo = 0)
  #capture.output({
    opt2 <- rstar(simu.des[[k]], thetainit = simu.pars.v, floglik = llik.fun,
                  fpsi = psi.fun,  psival = psi.fun(simu.pars.v),
                  datagen = gendat.fun, seed = 22, constr.opt = "solnp", R=250)
    opt1 <- optim(init, \(x) -opt.fn(x), method = "BFGS", hessian = TRUE)
  #}, type = c("message"))
  lik.vals[k] <- opt1$value
  lik.vals2[k] <- llik.fun(opt2$theta.hat, simu.des[[k]])
  par.stime[, k] <- crr.transform.par(opt1$par, inverse = TRUE)
  par.stime2[, k] <- crr.transform.par(opt2$par, inverse = TRUE)
  par.J[,,k] <- opt1$hessian
  #par.J2[,,k] <- opt2$hessian
  par.sd[, k] <- sqrt(diag(solve(opt1$hessian)))
  #par.sd2[, k] <- sqrt(diag(solve(opt2$hessian)))
}
dimnames(par.stime) <- dimnames(par.stime2) <- dimnames(par.sd) <- dimnames(par.sd2) <-
  list(pars = names(init), repl = seq_along(simu.des))
# fine simulazione
save.image(file.path(DIR, "sim1provv"))
