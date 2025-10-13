#' ---
#' title: studio di simulazione
#' subtitle: dati di Achana
#' format:
#'   pdf:
#'     fig-width: 12
#'     fig-height: 10
#' ---
rm(list = ls())
setwd("/home/marco/Nextcloud/tesi-ncrr/")
devtools::load_all()
library("likelihoodAsy")
library("tesi.ncrr")

CONFLVL <- .95
NSIM <- 250
VCOVTYPE <- "simple"

#DIR <- "../.." # per il markdown
DIR <- "/home/marco/output-tesi-morph/" # per l'esecuzione nel pacchetto
if (!dir.exists(DIR))
  dir.create(DIR)

#' # Simulazione basata sul problema di achana
#| warning: false
des <- ncrr.design(morphine)
des <- subset(des, which(sapply(des$design, \(d) (0 %in% d))))

simu.pars <- list(alpha = c(5.99639710290691, 1.51110953873727, -0.461698264404757), beta = c(0.551374273778573, 0.721813281778293, 0.751400418819617), mu0 = 37.4498851342351, sigma20 = 114.95716197389)

## simu.pars <- list(alpha = c(0.53118984013899, 1.0431973777787, 0.00434231242384523,
##                             2.36407165618289, 2.66293182986318, 2.7339581049579),
##                   beta = c(0.948918313002282,
##                            1.02425107553256, 1.06604309779969, 0.240340161234799,
##                            0.179383944934584, 0.179378353526596),
##                   mu0 = 0.81098898311333,
##                   sigma20 = 2.63212049308308,
##                   sigma2 = 5.69469982077026) # stime MV dai dati originali
simu.des <- simulate(des, nsim = NSIM, vcov.type = VCOVTYPE,
                     seed = 212, params = simu.pars)
simu.pars.v <- crr.join.par(simu.pars) |> crr.transform.par()

llik.fun <- get.llik.from.design(des, vcov.type = VCOVTYPE,
                                 echo = 0, use.data = TRUE)
gendat.fun <- function(data, theta) {
  theta <- crr.split.par(theta, length(data$treatments) - 1, transform = TRUE,
                         fixed = match.vcov.fixed(attr(data, "vcov.type")))
  suppressMessages(simulate.ncrr.design(data, params = theta, seed = 342)[[1]]) # estraggo solo il design, non tutta la lista
}
psi.fun <- function(theta) {
  theta[["sigma20"]]
}

save.image(file.path(DIR, "des.rda"))

#| eval: false
# simulazione
lik.vals <- lik.vals2 <- numeric(NSIM)
init <- getInitial(simu.des[[1]], vcov.type = VCOVTYPE)
par.J <- par.J2 <- array(NA, c(length(init), length(init), NSIM))
par.h0 <- par.stime <- par.stime2 <- par.sd <- par.sd2 <- matrix(NA, length(init), NSIM)
psi.r <- psi.rs <- psi.stime <- psi.sd <- numeric(NSIM)
for (k in seq_len(NSIM)) {
  message(sprintf("%.2f%%\r", k / NSIM * 100))
  # ------ OTTIMIZZAZIONE ----
  opt1 <- optim(init, \(x) -llik.fun(x, data = simu.des[[k]]),
                method = "BFGS", hessian = TRUE)

  for (R in as.integer(400*exp(1:2))){
    message("--> R = ", R)
    opt2 <- try(
      crr.rstar(simu.des[[k]], thetainit = init, floglik = llik.fun,
                fpsi = psi.fun,  psival = psi.fun(init), datagen = gendat.fun,
                seed = c(17980L, 31642L, 8590L, 104005L, 1543L, 257L),
                R = R, parallel = TRUE, trace = Inf)
    )
    if (!inherits(opt2, "try-error") && is.finite(opt2$rs)) break
  }

  saveRDS(list(optim = opt1, likasy = opt2),
          file = file.path(DIR, paste0("opt_", k, "_", as.integer(Sys.time()), ".rds")))
  # ------ REGISTRAZIONE RISULTATI ----
  lik.vals[k] <- opt1$value
  par.stime[, k] <- crr.transform.par(opt1$par, inverse = TRUE)
  par.J[,,k] <- opt1$hessian
  par.sd[, k] <- sqrt(diag(solve(opt1$hessian)))
  if (inherits(opt2, "try-error")) next
  lik.vals2[k] <- llik.fun(opt2$theta.hat, simu.des[[k]])
  par.stime2[, k] <- crr.transform.par(opt2$theta.hat, inverse = TRUE)
  par.J2[,,k] <- opt2$info.hat
  par.sd2[, k] <- opt2$se.theta.hat
  par.h0[, k] <- opt2$theta.hyp
  # stime di psi
  psi.r[k] <- opt2$r
  psi.rs[k] <- opt2$rs
  psi.stime[k] <- opt2$psi.hat
  psi.sd[k] <- opt2$se.psi.hat
}
dimnames(par.h0) <- dimnames(par.stime) <- dimnames(par.stime2) <-
  dimnames(par.sd) <- dimnames(par.sd2) <-
  list(pars = names(init), repl = seq_along(simu.des))
# fine simulazione
save.image(file.path(DIR, "sim1provv"))
