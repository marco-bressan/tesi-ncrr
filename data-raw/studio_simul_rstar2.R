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

#DIR <- "../.." # per il markdown
DIR <- "/home/marco/Nextcloud/output-tesi-morph/"
#DIR <- "../output-tesi" # per l'esecuzione nel pacchetto
load(file.path(DIR, "des.rda"))

if (file.exists(sim.file <- file.path(DIR, "sim1provv"))) {
  load(file.path)
} else {
  load(file.path(DIR, "des.rda"))
  #save(simu.des, file = file.path(DIR, "des.rda"))
  # risultati simulazione
  lik.vals <- lik.vals2 <- numeric(NSIM)
  init <- getInitial(simu.des[[1]], vcov.type = VCOVTYPE)
  par.J <- par.J2 <- array(NA, c(length(init), length(init), NSIM))
  par.h0 <- par.stime <- par.stime2 <- par.sd <- par.sd2 <- matrix(NA, length(init), NSIM)
  psi.r <- psi.rs <- psi.stime <- psi.sd <- numeric(NSIM)
  fs <- list.files(DIR, pattern = "^opt_")
  fnames <- do.call(rbind, strsplit(fs, "[._]"))
  dedup <- tapply(fnames[, 3], fnames[, 2], \(x) min(as.integer(x)))
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
    psi.r[k] <- obj$likasy$r
    psi.rs[k] <- obj$likasy$rs
    psi.stime[k] <- obj$likasy$psi.hat
    psi.sd[k] <- obj$likasy$se.psi.hat
    print(obj$likasy$rs)
  }
  dimnames(par.h0) <- dimnames(par.stime) <- dimnames(par.stime2) <-
    dimnames(par.sd) <- dimnames(par.sd2) <-
    list(pars = names(init), repl = seq_along(simu.des))
  save.image(file.path(DIR, "sim1provv"))
}


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


all.equal(par.stime[!is.na(par.stime)], par.stime2[!is.na(par.stime2)])

apply(par.stime2 - par.stime, 1, identity, simplify = FALSE) |>
  boxplot()
boxplot(lik.vals2 + lik.vals)

all.equal(par.J[-which(is.na(par.J2))], par.J2[-which(is.na(par.J2))])


# normalità di r e rstar

#par(mfrow = c(1, 3))
curve(dnorm(x), from = -10, to = 50, col = "red")
psi.r |> density() |> lines()
psi.rs |> na.omit() |> density() |> lines(col = "blue")



#| fig-cap: "Intervalli di confidenza Wald al 95% per i parametri (scala trasformata)"
#| fig-asp: 2
par(mfrow = c(4, 2))
for (j in seq_len(nrow(par.stime))) {
  print(pn <- names(simu.pars.v)[j])
  stima <- .ptrans(par.stime[j, ], pn, inverse = FALSE)
  ic <- cbind(stima,
              stima - qnorm(.975) * par.sd[j, ],
              stima + qnorm(.975) * par.sd[j, ])
  ic[is.infinite(ic)] <- sign(ic[is.infinite(ic)]) * 1000
  #if (pn == "sigma2") browser()
  par.vero.trans <- simu.pars.v[j]
  fuori <- par.vero.trans < ic[, 2] | par.vero.trans > ic[, 3]
  plot(ic[, 1], type = "n", main = pn,
       ylim = quantile(ic, c(.01, .99), na.rm = TRUE),
       sub = sprintf("Copertura empirica: %.2f%%",
                     mean(1 - fuori, na.rm = TRUE) * 100))
  segments(x0 = seq_len(nrow(ic)), y0 = ic[, 2], y1 = ic[, 3], lwd = 1.5,
           col = 1 + fuori)
  points(ic[, 1], pch = 16, col = 1 + fuori)
  abline(h = par.vero.trans, lty = 2, col = "blue")
}
