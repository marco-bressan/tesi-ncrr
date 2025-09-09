#' ---
#' title: studio di simulazione
#' subtitle: dati di Achana
#' format:
#'   pdf:
#'     fig-width: 12
#'     fig-height: 10
#' ---
rm(list = ls())
#setwd("/home/marco/Documenti/tesi-ncrr/")
devtools::load_all()
library("likelihoodAsy")
#library("tesi.ncrr")

CONFLVL <- .95
NSIM <- 250
VCOVTYPE <- "achana"

#DIR <- "/home/marco/output-tesi"
#DIR <- "../.." # per il markdown
DIR <- "../output-tesi/" # per l'esecuzione nel pacchetto

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
simu.pars.v <- tesi.ncrr:::crr.join.par(simu.pars) |> tesi.ncrr:::crr.transform.par()
gendat.fun <- function(data, theta) {
  theta <- tesi.ncrr:::crr.split.par(theta, length(data$treatments) - 1, transform = TRUE,
                         fixed = tesi.ncrr:::match.vcov.fixed(attr(data, "vcov.type")))
  suppressMessages(simulate(data, params = theta, seed = 342)[[1]]) # estraggo solo il design, non tutta la lista
}


psi.fun <- function(theta) {
  theta[["beta5"]]
}

rp.stat <- function(dati.gen, psi0, init, param = match("beta5", names(init)), ...) {
  # psi par d'interesse, lam di disturbo
  opt.theta <- nlminb(init, \(x) -llik.fun(x, dati.gen))
  theta.hat <- opt.theta$par
  psi.hat <- theta.hat[param]
  lam0 <- nlminb(theta.hat[-param], \(x) {
    z <- theta.hat
    z[param] <- psi0
    z[-param] <- x
    -llik.fun(z, dati.gen)
  })$par
  theta.hat[param] <- psi0
  theta.hat[-param] <- lam0
  lp0 <- llik.fun(theta.hat, dati.gen)
  if (-opt.theta$objective < lp0) browser()
  rp <- sign(psi.hat - psi0) * sqrt(2) * sqrt(-opt.theta$objective - lp0)
  structure(rp, theta = opt.theta)
}



#| eval: false
# simulazione
lik.vals <- lik.vals2 <- numeric(NSIM)
init <- getInitial(simu.des[[1]], vcov.type = VCOVTYPE)
par.J <- par.J2 <- array(NA, c(length(init), length(init), NSIM))
par.h0 <- par.stime <- par.stime2 <- par.sd <- par.sd2 <- matrix(NA, length(init), NSIM)
psi.rs <- psi.stime <- psi.sd <- numeric(NSIM)

for (k in 1) {
  k <- 1
  message(sprintf("%.2f%%\r", k / NSIM * 100))
  llik.fun <- get.llik.from.design(simu.des[[k]], vcov.type = VCOVTYPE, echo = 0)
  # ------ OTTIMIZZAZIONE ----
  opt1 <- crr.rstar(simu.des[[k]], thetainit = init, floglik = llik.fun,
                    fpsi = psi.fun,  psival = psi.fun(init),
                    datagen = gendat.fun, ronly = TRUE)
  unclass(opt1)
  boot.rp <- crr.boot(simu.des[[k]], rp.stat, R = 30, sim = "parametric",
                      ran.gen = gendat.fun, mle = simu.pars.v, parallel = FALSE,
                      seed = c(1998135100L, 2044097286L, 1091132551L, 966088075L, 1553350452L,
                               1303502678L),
                      psi0 = simu.pars.v[11], init = init, param = 11)
  boot.rp


  #saveRDS(list(optim = opt1, likasy = opt2),
  #        file = file.path(DIR, paste0("opt_", k, "_", as.integer(Sys.time()), ".rds")))
  # ------ REGISTRAZIONE RISULTATI ----
  ## lik.vals[k] <- opt1$value
  ## par.stime[, k] <- tesi.ncrr:::crr.transform.par(opt1$par, inverse = TRUE)
  ## par.J[,,k] <- opt1$hessian
  ## par.sd[, k] <- sqrt(diag(solve(opt1$hessian)))
  ## if (inherits(opt2, "try-error")) next
  ## lik.vals2[k] <- llik.fun(opt2$theta.hat, simu.des[[k]])
  ## par.stime2[, k] <- tesi.ncrr:::crr.transform.par(opt2$theta.hat, inverse = TRUE)
  ## par.J2[,,k] <- opt2$info.hat
  ## par.sd2[, k] <- opt2$se.theta.hat
  ## par.h0[, k] <- opt2$theta.hyp
  ## # stime di psi
  ## psi.rs[k] <- opt2$rs
  ## psi.stime[k] <- opt2$psi.hat
  ## psi.sd[k] <- opt2$se.psi.hat
}

dimnames(par.h0) <- dimnames(par.stime) <- dimnames(par.stime2) <-
  dimnames(par.sd) <- dimnames(par.sd2) <-
  list(pars = names(init), repl = seq_along(simu.des))
# fine simulazione
save.image(file.path(DIR, "sim1provv"))
