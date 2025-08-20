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
#DIR <- "../.." # per il markdown
DIR <- ".." # per l'esecuzione nel pacchetto

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
simu.des <- simulate(ncrr.design(smoke.alarm), vcov.type = VCOVTYPE,
                     nsim = NSIM, params = simu.pars)
#| eval: false
# simulazione
lik.vals <- lik.vals2 <- numeric(NSIM)
init <- getInitial(simu.des[[1]], vcov.type = VCOVTYPE)
par.J <- par.J2 <- array(NA, c(length(init), length(init), NSIM))
par.stime <- par.stime2 <- par.sd <- par.sd2 <- matrix(NA, length(init), NSIM)
for (k in seq_len(NSIM)) {
  cat(sprintf("%.2f%%\r", k / NSIM * 100))
  opt.fn <- get.llik.from.design(simu.des[[k]], vcov.type = VCOVTYPE, echo = 0)
  capture.output({
    #opt0 <- optim(init, \(x) -opt.fn(x), method = "Nelder-Mead", hessian = TRUE,
    #              control = list(maxit = 1000))
    #opt2 <- optim(opt0$par, \(x) -opt.fn(x), method = "BFGS", hessian = TRUE)
    opt1 <- optim(init, \(x) -opt.fn(x), method = "BFGS", hessian = TRUE)
  }, type = c("message"))
  lik.vals[k] <- opt1$value
  #lik.vals2[k] <- opt2$value
  par.stime[, k] <- crr.transform.par(opt1$par, inverse = TRUE)
  #par.stime2[, k] <- crr.transform.par(opt2$par, inverse = TRUE)
  par.J[,,k] <- opt1$hessian
  #par.J2[,,k] <- opt2$hessian
  par.sd[, k] <- sqrt(diag(solve(opt1$hessian)))
  #par.sd2[, k] <- sqrt(diag(solve(opt2$hessian)))
}
dimnames(par.stime) <- dimnames(par.stime2) <- dimnames(par.sd) <- dimnames(par.sd2) <-
  list(pars = names(init), repl = seq_along(simu.des))
# fine simulazione
save.image(file.path(DIR, "sim1provv"))

#| tbl-cap: "Confronto tra i parametri veri e quelli stimati con diversi
#| algoritmi di ottimizzazione"
load(file.path(DIR, "sim1.rda"))
load(file.path(DIR, "sim1provv"))
par.veri.v <- crr.join.par(simu.pars)
cbind(apply(par.stime, 1, median), apply(par.stime2, 1, median), par.veri.v) |>
  knitr::kable(col.names = c("BFGS", "Nelder-Mead + BFGS", "Par. veri"))

#| tbl-cap: "Confronto tra i valori ottimi della f.o. (-log(verosimiglianza))
#| con diversi algoritmi di ottimizzazione"
tab.lik <- cbind(1:250, lik.vals,lik.vals2,
                 round((lik.vals - lik.vals2)/lik.vals*100,2))
tab.lik <- cbind(head(tab.lik, 125), tail(tab.lik, 125))
knitr::kable(tab.lik, col.names = rep(c("# Repl.", "BFGS", "Nelder-Mead + BFGS", "Diff %"), 2))

#| fig-cap: >
#|   Boxplot delle stime di massima verosimiglianza dei parametri
#|   nelle varie replicazioni. Il "$+$" indica il valore dei veri parametri.
apply(par.stime, 1, identity, simplify = FALSE) |>
  boxplot()
points(seq_len(nrow(par.stime)), par.veri.v, cex = 2.5, col = "blue", pch = 3)


#| fig-cap: "Intervalli di confidenza Wald al 95% per i parametri (scala trasformata)"
#| fig-asp: 1.2
par(mfrow = c(4, 4))
for (j in seq_len(nrow(par.stime))) {
  print(pn <- names(par.veri.v)[j])
  stima <- .ptrans(par.stime[j, ], pn, inverse = FALSE)
  ic <- cbind(stima,
              stima - qnorm(.975) * par.sd[j, ],
              stima + qnorm(.975) * par.sd[j, ])
  ic[is.infinite(ic)] <- sign(ic[is.infinite(ic)]) * 1000
  #if (pn == "sigma2") browser()
  par.vero.trans <- .ptrans(par.veri.v[j], pn, FALSE)
  fuori <- par.vero.trans < ic[, 2] | par.vero.trans > ic[, 3]
  plot(ic[, 1], type = "n", main = pn,
       ylim = quantile(ic, c(.01, .99), na.rm = TRUE),
       sub = sprintf("Copertura empirica: %.2f%%",
                     mean(1 - fuori, na.rm = TRUE) * 100))
  segments(x0 = seq_len(nrow(ic)), y0 = ic[, 2], y1 = ic[, 3], lwd = 1,
           col = 1 + fuori)
  points(ic[, 1], pch = 16, col = 1 + fuori)
  abline(h = par.vero.trans, lty = 2, col = "blue")
}

#| fig-cap: "Intervalli di confidenza Wald al 95% per i parametri (scala originale)"
#| fig-asp: 1.5
par(mfrow = c(4, 4))
for (j in seq_len(nrow(par.stime))) {
  print(pn <- names(par.veri.v)[j])
  stima <- .ptrans(par.stime[j, ], pn, inverse = FALSE)
  ic <- cbind(stima,
              stima - qnorm(.975) * par.sd[j, ],
              stima + qnorm(.975) * par.sd[j, ])
  ic[] <- .ptrans(ic, pn, inverse = TRUE)
  ic[is.infinite(ic)] <- sign(ic[is.infinite(ic)]) * 1000
  #if (pn == "sigma2") browser()
  fuori <- par.veri.v[j] < ic[, 2] | par.veri.v[j] > ic[, 3]
  plot(ic[, 1], type = "n", main = pn,
       ylim = quantile(ic, c(.01, .99), na.rm = TRUE),
       sub = sprintf("Copertura empirica: %.2f%%",
                     mean(1 - fuori, na.rm = TRUE) * 100))
  segments(x0 = seq_len(nrow(ic)), y0 = ic[, 2], y1 = ic[, 3], lwd = 1,
           col = 1 + fuori)
  points(ic[, 1], pch = 16, col = 1 + fuori)
  abline(h = par.veri.v[j], lty = 2, col = "blue")
}

#' # intervalli LRT
#'
llik <- get.llik.from.design(simu.des[[1]], vcov.type = VCOVTYPE)
likfun <- function(x, i, j, y = par.stime[, i]) {
  y <- crr.transform.par(y) # da parametri originali a trasformati
  y[j] <- x
  if (any(is.nan(y))) stop("Nan (llik)")#browser()
  llik(y, crr.get.theta(simu.des[[i]], raw = TRUE),
       crr.get.Gamma(simu.des[[i]], raw = TRUE))
}
trvfun <- function(x, i, j, y = par.stime[, i], ly = lik.vals[i], J = par.J[, , i]) {
  logL <- \(z) llik(z, crr.get.theta(simu.des[[i]], raw = TRUE),
                   crr.get.Gamma(simu.des[[i]], raw = TRUE))
  #browser()
  y.old <- y <- crr.transform.par(y) # da parametri originali a trasformati
  y[-j] <- y[-j] + c(solve(J[-j, -j]) %*% J[-j, j] * (y[j] - x)) # parametro di disturbo
  #if (abs(x - y[j]) < 0.1) browser()
  y[j] <- x # parametro d'interesse
  if (any(!is.finite(y))) browser()
  ly + logL(y) + qchisq(CONFLVL, 1) / 2
}

#| eval: false
# calcolo degli intervalli TRV
trv <- array(NA, c(dim(par.stime), 2))
cc <- 0
dimnames(trv) <- append(dimnames(par.stime), list(side = c("lower", "upper")))
for (i in 1:ncol(par.stime)) {
  llik <- get.llik.from.design(simu.des[[i]], vcov.type = VCOVTYPE)
  par.trans <- crr.transform.par(par.stime[, i])
  for (j in 1:nrow(par.stime)) {
    pn <- rownames(par.stime)[j]
    cat(sprintf("%.2f%%\r", (cc <- cc + 1) / prod(dim(par.stime)) * 100))
    cerca.in <- par.trans[j] + c(-10, 10) *
      ifelse(is.finite(par.sd[j, i]), par.sd[j, i], par.stime[j, i])
    #max.nuovo <- optimise(\(x) trvfun(x, i, j), cerca.in, maximum = TRUE)$maximum
    suppressMessages({
      trvci <- try({
        c(uniroot(\(x) trvfun(x, i, j), c(cerca.in[1], par.trans[j]))$root,
          uniroot(\(x) trvfun(x, i, j), c(par.trans[j], cerca.in[2]))$root)
      }, silent = TRUE)
    })#, type = c("message"))
    if (!inherits(trvci, "try-error")) {
      trv[j, i, ] <- .ptrans(trvci, pn, TRUE)
    } else {
      browser()
      plot(Vectorize(\(x) trvfun(x, i, j)), from = cerca.in[1], to = cerca.in[2])
      abline(v = cerca.in, col="gray30", lty = 2)
      abline(h=0, col="red")
      abline(v=par.trans[j], col = "blue")
    }
  }
}
save(trv, file = file.path(DIR, "sim1_trv.rda"))

#| eval: true
load(file.path(DIR, "sim1_trv.rda"))
apply(trv, c(1, 3), \(x) mean(is.finite(x)))


#| fig-asp: 0.9
par.veri.trans <- crr.transform.par(par.veri.v)
set.seed(12)
par(mfrow = c(4, 4))
for (i in sample(seq_along(simu.des), 4)) {
  llik <- get.llik.from.design(simu.des[[i]], vcov.type = VCOVTYPE)
  par.trans <- crr.transform.par(par.stime[, i])
  for (j in sample(1:nrow(par.stime), 4)) {
    pn <- rownames(par.stime)[j]
    cerca.in <- par.trans[j] + c(-5, 5) *
      ifelse(is.finite(par.sd[j, i]), par.sd[j, i], par.stime[j, i])
    #max.nuovo <- optimise(\(x) trvfun(x, i, j), cerca.in, maximum = TRUE)$maximum
    curve(Vectorize(\(z) trvfun(z, i, j))(x), from = cerca.in[1], to = cerca.in[2],
          main = bquote(W[paste("repl.", .(i))](.(pn))))
    abline(v = c(par.veri.trans[j], par.trans[j], trv[j, i, ]),
           col = c("green", "blue", "red", "red"), lty = c(1, 1, 2, 2))
    abline(h = 0, col = "gray30", lty = 2)
  }
}

#| fig-cap: "Intervalli di confidenza al 95% basati sulla verosimiglianza per i parametri (scala originale)."
#| fig-asp: 2
par(mfrow = c(5, 3))
for (j in seq_len(nrow(par.stime))) {
  ic <- cbind(par.stime[j, ], trv[j, , ])
  #if (rownames(par.stime)[j] == "mu0") browser()
  fuori <- par.veri.v[j] < ic[, 2] | par.veri.v[j] > ic[, 3]
  plot(ic[, 1], type = "n", main = rownames(par.stime)[j],
       ylim = range(ic, na.rm = TRUE),
       sub = sprintf("Copertura empirica: %.2f%%",
                     mean(1 - fuori, na.rm = TRUE) * 100))
  segments(x0 = seq_len(nrow(ic)), y0 = ic[, 2], y1 = ic[, 3], lwd = 1,
           col = 1 + fuori)
  points(ic[, 1], pch = 16, col = 1 + fuori)
  abline(h = par.veri.v[j], lty = 2, col = "blue")
