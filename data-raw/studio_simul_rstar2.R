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

CONFLVL <- .95
NSIM <- 250
VCOVTYPE <- "achana"

DIR <- r"{C:\Users\marco\Nextcloud\output-tesi}"
#DIR <- "../.." # per il markdown
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
simu.pars.v <- tesi.ncrr:::crr.join.par(simu.pars) |> tesi.ncrr:::crr.transform.par()
llik.fun <- function(theta, data) {
  y = tesi.ncrr:::crr.get.theta(data, raw = TRUE)
  Gamma = tesi.ncrr:::crr.get.Gamma(data, raw = TRUE)
  L <- llik(theta, y = y, Gamma = Gamma)
  # if (!is.finite(L)) {
  #   environment(llik)$echo <- 3
  #   browser()
  #   llik(theta, y = y, Gamma = Gamma)
  # }
  return(L)
}
gendat.fun <- function(theta, data) {
  theta <- tesi.ncrr:::crr.split.par(theta, length(data$treatments) - 1, transform = TRUE,
                         fixed = tesi.ncrr:::match.vcov.fixed(attr(data, "vcov.type")))
  suppressMessages(simulate(data, params = theta, seed = 342)[[1]]) # estraggo solo il design, non tutta la lista
}
psi.fun <- function(theta) {
  theta[["beta5"]]
}

save(simu.des, file = file.path(DIR, "des.rda"))

#| eval: false
# simulazione
lik.vals <- lik.vals2 <- numeric(NSIM)
init <- getInitial(simu.des[[1]], vcov.type = VCOVTYPE)
par.J <- par.J2 <- array(NA, c(length(init), length(init), NSIM))
par.h0 <- par.stime <- par.stime2 <- par.sd <- par.sd2 <- matrix(NA, length(init), NSIM)
psi.rs <- psi.stime <- psi.sd <- numeric(NSIM)

fs <- list.files(DIR, pattern = "^opt_")
fnames <- do.call(rbind, strsplit(fs, "[._]"))
dedup <- tapply(fnames[, 3], fnames[, 2], \(x) max(as.integer(x)))
dedup.i <- match(dedup, as.integer(fnames[, 3]))
fs <- fs[dedup.i]
fnames <- fnames[dedup.i, ]
for (i in seq_along(fs)) {
  # --- RECUPERO DA DISCO ---
  k <- as.integer(fnames[i, 2])
  llik <- get.llik.from.design(simu.des[[k]], vcov.type = VCOVTYPE, echo = 0)
  obj <- readRDS(file.path(DIR, fs[i]))
  # ------ REGISTRAZIONE RISULTATI ----
  lik.vals[k] <- obj$optim$value
  par.stime[, k] <- tesi.ncrr:::crr.transform.par(obj$optim$par, inverse = TRUE)
  par.J[,,k] <- obj$optim$hessian
  par.sd[, k] <- sqrt(diag(solve(obj$optim$hessian)))
  if (inherits(obj$likasy, "try-error")) next
  lik.vals2[k] <- llik.fun(obj$likasy$theta.hat, simu.des[[k]])
  par.stime2[, k] <- tesi.ncrr:::crr.transform.par(obj$likasy$theta.hat, inverse = TRUE)
  par.J2[,,k] <- obj$likasy$info.hat
  par.sd2[, k] <- obj$likasy$se.theta.hat
  par.h0[, k] <- obj$likasy$theta.hyp
  # stime di psi
  psi.rs[k] <- obj$likasy$rs
  psi.stime[k] <- obj$likasy$psi.hat
  psi.sd[k] <- obj$likasy$se.psi.hat
}

dimnames(par.h0) <- dimnames(par.stime) <- dimnames(par.stime2) <-
  dimnames(par.sd) <- dimnames(par.sd2) <-
  list(pars = names(init), repl = seq_along(simu.des))

save.image(file.path(DIR, "sim1provv"))

par.stime


# confronti grafici

#| fig-cap: >
#|   Boxplot delle stime di massima verosimiglianza dei parametri
#|   nelle varie replicazioni. Il "$+$" indica il valore dei veri parametri.
par(mfrow = c(1, 2))
apply(par.stime, 1, identity, simplify = FALSE) |>
  boxplot()
points(seq_len(nrow(par.stime)), do.call(c, simu.pars), cex = 2.5, col = "blue", pch = 4)
apply(par.stime2, 1, identity, simplify = FALSE) |>
  boxplot()
points(seq_len(nrow(par.stime)), do.call(c, simu.pars), cex = 2.5, col = "blue", pch = 4)

apply(par.stime2 - par.stime, 1, identity, simplify = FALSE) |>
  boxplot()
cbind(lik.vals2, -lik.vals)
lik.vals2 + lik.vals |> boxplot()

all.equal(par.J[-which(is.na(par.J2))], par.J2[-which(is.na(par.J2))])
