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

smoke.alarm

smoke.alarm$tik <- with(smoke.alarm, log(rik) - log(nik - rik))
names(smoke.alarm)

# specifico il design della meta-analisi
des <- ncrr.design(smoke.alarm)

#' Nota: la funzione avverte che i baseline != 0 non sono implementati.
#' In realtà, per `vcov="achana"` lo sarebbero ma per confronto con gli
#' altri metodi al momento decido di toglierli dal dataset

# solo i design che contengono lo zero
des <- subset(des, which(sapply(des$design, \(x) 0 %in% x)))

opt.fn <- get.llik.from.design(des, vcov.type = "achana")
opt1 <- optim(getInitial(des, vcov.type = "achana"),
              \(x) -opt.fn(x), method = "BFGS")
crr.split.par(opt1$par, 5, transform = TRUE)

opt2 <- optim(opt1$par, \(x) -opt.fn(x), method = "Nelder-Mead")
crr.split.par(opt2$par, 5, transform = TRUE)

#'
#' Parametri (I colonna) con relativi s.e.
#'
cbind(opt2$par,
      optimHess(opt2$par, \(x) -opt.fn(x)) |>
        solve() |>
        diag() |>
        sqrt()
      )

opt.fn <- get.llik.from.design(des, vcov.type = "normal", echo = 3)
as.list(environment(opt.fn))
opt1 <- optim(getInitial(des, vcov.type = "normal"), \(x) -opt.fn(x), method = "BFGS")
crr.split.par(opt1$par, 5, transform = TRUE)

opt2 <- optim(opt1$par, \(x) -opt.fn(x), method = "Nelder-Mead")
crr.split.par(opt2$par, 5, transform = TRUE)

# escludo rho, che non viene stimato in quanto ci sono solo studi
# "bivariati"
opt2h <- optimHess(opt1$par, \(x) -opt.fn(x))
cbind(opt1$par[-13],
      opt2h[-13, ][, -13] |>
        solve() |>
        diag() |>
        sqrt()
      )

opt.fn <- get.llik.from.design(des, vcov.type = "simplified")
as.list(environment(opt.fn))
opt1 <- optim(getInitial(des, vcov.type = "simplified"),
              \(x) -opt.fn(x), method = "BFGS")
crr.split.par(opt1$par, 5, transform = TRUE)

opt2h <- optimHess(opt1$par, \(x) -opt.fn(x))
cbind(opt1$par[-13],
      opt2h[-13, ][, -13] |>
        solve() |>
        diag() |>
        sqrt()
      )



#' # Esempio 2
#'
#' In questo caso gli *outcome* della NMA sono numerici, quindi il
#' modello può essere considerato normale e le varianze associate
#' a ciascuno studio (riportate esplicitamente nello studio) sono
#' inserite direttamente nella matrice $\hat \Gamma = diag(s^2_{ik})$.

# specifico il design della meta-analisi
des2 <- ncrr.design(morphine)
str(des2)
des2 <- subset(des2, which(sapply(des$design, \(x) 0 %in% x)))

opt.fn <- get.llik.from.design(des2, vcov.type = "achana")
opt1 <- optim(getInitial(des2, vcov.type = "achana"),
              \(x) -opt.fn(x), method = "BFGS")
crr.split.par(opt1$par, 3, transform = TRUE)
cbind(opt1$par,
      optimHess(opt1$par, \(x) -opt.fn(x)) |>
        solve() |>
        diag() |>
        sqrt()
      )

opt2.fn <- get.llik.from.design(des2, vcov.type = "normal")
opt21 <- optim(getInitial(des2, vcov.type = "normal"),
              \(x) -opt2.fn(x), method = "BFGS")
crr.split.par(opt21$par, 3, transform = TRUE)
cbind(opt21$par[-9],
      optimHess(opt21$par, \(x) -opt2.fn(x))[-9, ][, -9] |>
        solve() |>
        diag() |>
        sqrt()
      )


opt.fn <- get.llik.from.design(des, vcov.type = "normal")
as.list(environment(opt.fn))
opt1 <- optim(getInitial(des, vcov.type = "normal"), \(x) -opt.fn(x), method = "BFGS")
crr.split.par(opt1$par, 5, transform = TRUE)

opt2 <- optim(opt1$par, \(x) -opt.fn(x), method = "Nelder-Mead")
crr.split.par(opt2$par, 5, transform = TRUE)

opt2h <- optimHess(opt1$par, \(x) -opt.fn(x))
cbind(opt1$par[-13],
      opt2h[-13, ][, -13] |>
        solve() |>
        diag() |>
        sqrt()
      )
