#' ---
#' title: "prova dell'ottimizzatore per la network crr"
#' subtitle: "Replicazione dell'esperimento di Achana et al."
#' format:
#'   pdf:
#'     fig-width: 8
#'     fig-height: 6
#' execute:
#'   eval: true
#'   echo: false
#'   output: true
#' ---
#'
#' # Importazione dei dati
#'
#' Dati sull'impatto dell'educazione domiciliare per l'installazione di
#' un impianto di rilevazione del fumo.
rm(list = ls())
#| output: false
devtools::load_all(".")
data("smoke.alarm")
smoke.alarm$tik <- with(smoke.alarm, log(rik) - log(pmax(nik - rik, .001)))

# specifico il design della meta-analisi
des <- ncrr.design(smoke.alarm)

ini1 <- getInitial(des, vcov.type = "achana", rep = 1)
curr <- crr.split.par(ini1, np = 6, transform = TRUE, fixed = match.vcov.fixed("achana"))
Gamma <- crr.get.Gamma(des, raw = TRUE)
Y <- crr.get.theta(des, raw = TRUE)
n <- length(des$design)

# E-step
#TODO: non va bene per la network crr! deve essere generalizzato, ma ciò è possibile solo conoscendo
#esattamente come sono state calcolate tutte le quantità in gioco
curr.mu <- crr.get.mu(des, curr, raw = TRUE)
curr.Sigma <- crr.get.sigma(des, set.vcov.params(curr, vcov.type = "achana"), raw = TRUE)
Lambda.star <- mapply(\(Gi, Si) solve(solve(Gi) + solve(Si)), Gamma, curr.Sigma, SIMPLIFY = FALSE)
theta.star <- mapply(\(Gi, Si, Ssi, yi, mui) Ssi %*% (solve(Si) %*% mui + solve(Gi) %*% yi),
                     Gamma, curr.Sigma, Lambda.star, Y, curr.mu, SIMPLIFY = FALSE)
theta.m <- 1/n * Reduce("+", theta.star)
theta.S <- 1/n * Reduce("+", Lambda.star[-1], init = Lambda.star[[1]]) +
  1/n * Reduce("+", lapply(theta.star, \(x) tcrossprod(x - theta.m)), c(0, 0))

# M-step



## Confronto con ML

#' Nota: la funzione avverte che i baseline != 0 non sono implementati.
#' In realtà, per `vcov="achana"` lo sarebbero ma per confronto con gli
#' altri metodi al momento decido di toglierli dal dataset

# solo i design che contengono lo zero
#des <- subset(des, which(sapply(des$design, \(x) 0 %in% x)))

opt.fn <- get.llik.from.design(des, vcov.type = "achana", echo = 0)
opt1 <- optim(ini1,
              \(x) -opt.fn(x), method = "BFGS")
crr.split.par(opt1$par, 6, transform = TRUE, fixed = match.vcov.fixed("achana"))

#mettere come primo passo
opt2 <- optim(opt1$par, \(x) -opt.fn(x), method = "Nelder-Mead")
crr.split.par(opt2$par, 6, transform = TRUE, fixed = match.vcov.fixed("achana"))

#'
#' Parametri (I colonna) con relativi s.e.
#'
cbind(par = crr.transform.par(opt2$par, np = 6, split = FALSE,
                        fixed = match.vcov.fixed("achana"), inverse = TRUE),
      stderr = optimHess(opt1$par, \(x) -opt.fn(x)) |>
        solve() |>
        diag() |>
        sqrt()
      ) |>
  round(4)
