#' ---
#' title: studio di simulazione
#' subtitle: dati di Achana
#' format:
#'   pdf:
#'     fig-width: 12
#'     fig-height: 10
#' ---
rm(list = ls())
#setwd("/home/marco/UNIPD Scienze Statistiche/tesi-ncrr/")
devtools::load_all()

CONFLVL <- .95
NSIM <- 250
VCOVTYPE <- "achana"
#DIR <- "../.." # per il markdown
DIR <- ".." # per l'esecuzione nel pacchetto

#' # Simulazione basata sul problema di achana
#| warning: false
simu.pars <- list(alpha = rep(-1:1, length = 6), beta = rep(1:3, length = 6),
                  mu0 = 1, sigma20 = 0.5, sigma2 = 1)
simu.des <- do.call(simulate,
                    append(list(ncrr.design(smoke.alarm), vcov.type = VCOVTYPE,
                                nsim = NSIM),
                           simu.pars))

#| eval: false
# simulazione
lik.vals <- lik.vals2 <- numeric(NSIM)
init <- getInitial(simu.des[[1]], vcov.type = VCOVTYPE)
par.stime <- par.stime2 <- par.sd <- par.sd2 <- matrix(NA, length(init), NSIM)
for (k in seq_len(NSIM)) {
  cat(sprintf("%.2f%%\r", k / NSIM * 100))
  opt.fn <- get.llik.from.design(simu.des[[k]], vcov.type = VCOVTYPE, echo = 0)
  capture.output({
    opt0 <- optim(init, \(x) -opt.fn(x), method = "Nelder-Mead", hessian = TRUE,
                  control = list(maxit = 1000))
    opt2 <- optim(opt0$par, \(x) -opt.fn(x), method = "BFGS", hessian = TRUE)
    opt1 <- optim(init, \(x) -opt.fn(x), method = "BFGS", hessian = TRUE)
  }, type = c("message"))
  lik.vals[k] <- opt1$value
  lik.vals2[k] <- opt2$value
  par.stime[, k] <- crr.transform.par(opt1$par, inverse = TRUE)
  par.stime2[, k] <- crr.transform.par(opt2$par, inverse = TRUE)
  par.sd[, k] <- sqrt(diag(solve(opt1$hessian)))
  par.sd2[, k] <- sqrt(diag(solve(opt2$hessian)))
}
dimnames(par.stime) <- dimnames(par.stime2) <- dimnames(par.sd) <- dimnames(par.sd2) <-
  list(pars = names(init), repl = seq_along(simu.des))
# fine simulazione

#' # Intervalli alla Wald

#| tbl-cap: "Confronto tra i parametri veri e quelli stimati con diversi
#| algoritmi di ottimizzazione"
load(file.path(DIR, "sim12.rda"))
par.veri.v <- crr.join.par(simu.pars)
cbind(apply(par.stime, 1, median), apply(par.stime2, 1, median), par.veri.v) |>
  knitr::kable(col.names = c("BFGS", "Nelder-Mead + BFGS", "Par. veri"))

#| tbl-cap: "Confronto tra i valori ottimi della f.o. (-log(verosimiglianza))
#| con diversi algoritmi di ottimizzazione"
tab.lik <- cbind(1:250, lik.vals,lik.vals2,
                 round((lik.vals - lik.vals2)/lik.vals*100,2))
tab.lik <- cbind(head(tab.lik, 125), tail(tab.lik, 125))
knitr::kable(tab.lik, col.names = rep(c("# Repl.", "BFGS", "Nelder-Mead + BFGS", "Diff %"), 2))

##na.cols <- col(par.sd)[c(which(!is.finite(par.sd2)), which(!is.finite(par.stime2)))]
## if (length(na.cols) > 0) {
##   par.stime2 <- par.stime2[, -na.cols]
##   par.sd2 <- par.sd2[, -na.cols]
## }



#| fig-cap: "Intervalli di confidenza Wald al 95% per i parametri (scala trasformata)"
#| fig-asp: 2
par(mfrow = c(5, 3))
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
  segments(x0 = seq_len(nrow(ic)), y0 = ic[, 2], y1 = ic[, 3], lwd = 1.5,
           col = 1 + fuori)
  points(ic[, 1], pch = 16, col = 1 + fuori)
  abline(h = par.vero.trans, lty = 2, col = "blue")
}

#| fig-cap: "Intervalli di confidenza Wald al 95% per i parametri (scala originale)"
#| fig-asp: 2
par(mfrow = c(5, 3))
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
  segments(x0 = seq_len(nrow(ic)), y0 = ic[, 2], y1 = ic[, 3], lwd = 1.5,
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
trvfun <- function(x, i, j, y = par.stime[, i], ly = lik.vals[i]) {
  y <- crr.transform.par(y) # da parametri originali a trasformati
  y[j] <- x
  if (any(is.nan(y))) stop("Nan (trv)")#browser()
  #browser()
  ly + llik(y, crr.get.theta(simu.des[[i]], raw = TRUE),
            crr.get.Gamma(simu.des[[i]], raw = TRUE)) +
    qchisq(CONFLVL, 1) / 2
}

#| eval: false
cc <- 0
errori.stima.mv <- array(NA, dim(par.stime))
par.stime22 <- par.stime
par.sd22 <- par.sd
for (i in 1:ncol(par.stime)) {
  par.trans <- crr.transform.par(par.stime[, i])
  for (j in 1:nrow(par.stime)) {
    pn <- rownames(par.stime)[j]
    cat(sprintf("%.2f%%\r", (cc <- cc + 1) / prod(dim(par.stime)) * 100))
    # NOTA: ottimizzazione univariata (scala trasformata) fissando gli altri parametri
    # è il modo SBAGLIATO di calcolare la verosimiglianza marginale!
    #if (j >= 14) browser()
    cerca.in <- par.trans[j] + c(-10, 10) *
      ifelse(is.finite(par.sd[j, i]), par.sd[j, i], par.stime[j, i])
    max.nuovo <- optimise(\(x) likfun(x, i, j), cerca.in, maximum = TRUE)$maximum
    # errori % commessi sulla scala originale
    errori.stima.mv[j, i] <- abs(par.trans[j] - max.nuovo) / abs(par.trans[j])
    par.trans[j] <- max.nuovo
  }
  par.stime22[, i] <- crr.transform.par(par.trans, inverse = TRUE)
  # ricalcolo hessiana
  valida <- par.trans > -100 & par.trans < 100
  par.sd22[valida, i] <- suppressWarnings(sqrt(diag(solve(
    optimHess(crr.transform.par(par.stime22[, i]),
              \(x) -llik(x, crr.get.theta(simu.des[[i]], raw = TRUE),
                         crr.get.Gamma(simu.des[[i]], raw = TRUE)))[valida, valida]
  ))))
  par.sd22[!valida, i] <- NA
  if (any(!valida)) {
    warning("Rilevata stima al bordo dello spazio parametrico per ",
            paste(names(par.trans)[!valida], collapse = ", "),
            sprintf(" (alla replicazione %i)", i))
  }
}
save(lik.vals,lik.vals2, par.stime, par.stime2, par.sd, par.sd2, par.sd22, par.stime22,
     errori.stima.mv, file = file.path(DIR, "sim12.rda"))


#| eval: true
#| fig-cap: >
#|   Campione casuale di 12 verosimiglianze marginali in diverse replicazioni.
#|   La linea continua identifica la verosimiglianza profilo fissando gli altri
#|   parametri in base alla stima di massima verosimiglianza ottenuta tramite
#|   Nelder-Mead + BFGS. La linea tratteggiata, invece, è la verosimiglianza profilo
#|   fissando i parametri ai loro valori ottimi dopo aver trovato le stime
#|   univariate. La linea blu indica la stima di massima verosimiglianza non
#|   vincolata; la linea verde indica il massimo vincolato. La scala usata è
#|   quella trasformata, ovvero: nessuna trasformazione per i parametri di
#|   posizione $\alpha_*, \beta_*, \mu_0$; trasformazione logaritmica per le
#|   varianze $\sigma_*$; trasformata di Fisher per la correlazione $\rho$.
#| fig-asp: 1.2
par(mfrow = c(3, 4))
# campionamento casuale delle stime univariate ottenute
# (mostro grafico solo in presenza di errori percentuali vistosi)
set.seed(3502)
idx <- sample(seq_along(errori.stima.mv), 12)
nomi <- rownames(par.stime22)[1 + idx %% nrow(par.stime22)]
replicazioni <- 1 + idx %/% nrow(par.stime22)
for (i in seq_along(idx)) {
  j <- match(nomi[i], rownames(par.stime22))
  k <- replicazioni[i]
  par.trans <- crr.transform.par(par.stime[, k])
  xrng <- par.trans[j] + c(-10, 10) *
    ifelse(is.finite(par.sd[j, k]), par.sd[j, k], par.stime[j, k])
  varsym <- as.name(nomi[i])
  # scelgo intervallo abbastanza ampio per essere sicuro
  xx <- seq(xrng[1], xrng[2], length.out = 100)
  likvals <- sapply(xx, \(x) likfun(x, k, j))
  likvals2 <- sapply(xx, \(x) likfun(x, k, j, y = par.stime22[, k]))
  # il grafico lo faccio sulla scala trasformata
  plot(xx, likvals, type = "l",
       xlab = paste("Parametro", if (grepl("sigma", nomi[i])) "(scala trasformata)"),
       main = bquote(log ~ L["repl." ~ .(k)](.(varsym))),
       ylim = quantile(c(likvals, likvals2), c(.05, 1)))
  lines(xx, likvals2, lty = 2)
  abline(v = par.trans[j], col = "blue")
  abline(v = .ptrans(par.stime22[j, k], nomi[i], FALSE), col = "green")
}

lik.vals22 <- numeric(length(lik.vals2))
for (i in seq_along(lik.vals22)) {
  lik.vals22[i] <- -llik(crr.transform.par(par.stime22[, i]),
                         crr.get.theta(simu.des[[i]], raw = TRUE),
                         crr.get.Gamma(simu.des[[i]], raw = TRUE))
}

cbind(
  apply(par.sd, 1, \(x) mean(is.finite(x))),
  apply(par.sd22, 1, \(x) mean(is.finite(x)))
)

#| fig-cap: "Confronto tra gli algoritmi di ottimizzazione"
boxplot(list("univariata" = lik.vals22,
             "BFGS" = lik.vals,
             "N-M + BFGS" = lik.vals22))

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
    }
  }
}
save(trv, file = file.path(DIR, "sim12_trv.rda"))

#| eval: true
load(file.path(DIR, "sim12_trv.rda"))
apply(trv, c(1, 3), \(x) mean(is.finite(x)))

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
  segments(x0 = seq_len(nrow(ic)), y0 = ic[, 2], y1 = ic[, 3], lwd = 1.5,
           col = 1 + fuori)
  points(ic[, 1], pch = 16, col = 1 + fuori)
  abline(h = par.veri.v[j], lty = 2, col = "blue")
}
