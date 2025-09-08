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

# specifico il design della meta-analisi
des <- ncrr.design(smoke.alarm)
#str(des)

#' Nota: la funzione avverte che i baseline != 0 non sono implementati.
#' In realtà, per `vcov="achana"` lo sarebbero ma per confronto con gli
#' altri metodi al momento decido di toglierli dal dataset

# solo i design che contengono lo zero
#des <- subset(des, which(sapply(des$design, \(x) 0 %in% x)))

confronta.gr <- function(x, fn, perf = FALSE) {
  numgr <- numDeriv::grad(fn, x)
  gr <- attr(fn, "score")(x)
  if (perf)
    perf <- microbenchmark::microbenchmark(
      numgr = numDeriv::grad(fn, x),
      gr = attr(fn, "score")(x),
      times = 20
    )
  print(perf)
  cat("\n========================\n")
  cat("Differenze:", paste("\t-", all.equal(numgr, gr)), sep = "\n")
  cat("\n========================\n")
  cbind(x = x, numgr = numgr, gr = gr)
}

opt.fn <- get.llik.from.design2(des, vcov.type = "achana", echo = 0)
opt1 <- optim(ini1 <- getInitial(des, vcov.type = "achana"),
              \(x) -opt.fn(x), method = "BFGS")
crr.split.par(opt1$par, 6, transform = TRUE, fixed = match.vcov.fixed("achana"))

#mettere come primo passo
ini12 <- getInitial(des, vcov.type = "achana", seed = 2)
confronta.gr(ini1, opt.fn, perf = F)
confronta.gr(ini12, opt.fn, perf = F)
confronta.gr(opt1$par, opt.fn, perf = F)

opt2 <- optim(ini1, \(x) -opt.fn(x), method = "BFGS", gr = \(x) -1e-7 * attr(opt.fn, "score")(x))
confronta.gr(opt2$par, opt.fn, perf = F)
nlminb(ini1, \(x) -opt.fn(x), gradient = \(x) -1e-7 * attr(opt.fn, "score")(x)) # false convergence
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

opt.fn <- get.llik.from.design(des, vcov.type = "normal", echo = 4, transform = TRUE)
opt1 <- optim(ini2 <- getInitial(des, vcov.type = "normal", transform = TRUE),
              \(x) -opt.fn(x), method = "BFGS")
#crr.split.par(opt1$par, 5, transform = TRUE)

opt2 <- optim(ini2, \(x) -opt.fn(x), method = "Nelder-Mead")
crr.split.par(opt2$par, 6, transform = TRUE)

opt2h <- optimHess(opt2$par, \(x) -opt.fn(x))
cbind(crr.transform.par(opt2$par, np = 6, split = FALSE, inverse = TRUE),
      opt2h |>
        solve() |>
        diag() |>
        sqrt()) |>
  round(4)

opt.fn <- get.llik.from.design(des, vcov.type = "simple")
opt1 <- optim(ini3 <- getInitial(des, vcov.type = "simple"),
              \(x) -opt.fn(x), method = "BFGS")
crr.split.par(opt1$par, 6, transform = TRUE, fixed = match.vcov.fixed("simple"))

opt2h <- optimHess(opt1$par, \(x) -opt.fn(x))
cbind(crr.transform.par(opt1$par, np = 6, split = FALSE,
                        inverse = TRUE, fixed = match.vcov.fixed("simple"))[-13],
      opt2h |>
        solve() |>
        diag() |>
        sqrt())






#' # Esempio 2
#'
#' In questo caso gli *outcome* della NMA sono numerici, quindi il
#' modello può essere considerato normale e le varianze associate
#' a ciascuno studio (riportate esplicitamente nello studio) sono
#' inserite direttamente nella matrice $\hat \Gamma = diag(s^2_{ik})$.

# specifico il design della meta-analisi
des2 <- ncrr.design(morphine)
#str(des2)
des2 <- subset(des2, which(sapply(des2$design, \(x) 0 %in% x)))

opt.fn <- get.llik.from.design(des2, vcov.type = "achana")
opt1 <- optim(getInitial(des2, vcov.type = "achana"),
              \(x) -opt.fn(x), method = "BFGS")
crr.split.par(opt1$par, 3, transform = TRUE, fixed = match.vcov.fixed("achana"))
cbind(crr.transform.par(opt1$par, inverse = TRUE, split = FALSE,
                        np = 3, fixed = match.vcov.fixed("achana")),
      optimHess(opt1$par, \(x) -opt.fn(x)) |>
        solve() |>
        diag() |>
        sqrt()
      ) |> round(4)

opt2.fn <- get.llik.from.design(des2, vcov.type = "normal")
opt21 <- optim(getInitial(des2, vcov.type = "normal"),
              \(x) -opt2.fn(x), method = "BFGS")
crr.split.par(opt21$par, 3, transform = TRUE)
cbind(crr.transform.par(opt21$par, inverse = TRUE, split = FALSE, np = 3),
      optimHess(opt21$par, \(x) -opt2.fn(x)) |>
        solve() |>
        diag() |>
        sqrt()) |> round(4)


opt.fn <- get.llik.from.design(des2, vcov.type = "equivar")
opt1 <- optim(getInitial(des2, vcov.type = "equivar"), \(x) -opt.fn(x), method = "BFGS")
crr.split.par(opt1$par, 3, transform = TRUE, fixed = match.vcov.fixed("equivar"))
opt2 <- optim(opt1$par, \(x) -opt.fn(x), method = "Nelder-Mead")
crr.split.par(opt2$par, 3, transform = TRUE, fixed = match.vcov.fixed("equivar"))

opt2h <- optimHess(opt1$par, \(x) -opt.fn(x))
cbind(crr.transform.par(opt2$par, np = 3, inverse = TRUE, split = FALSE,
                    fixed = match.vcov.fixed("equivar")),
      opt2h[-13, ][, -13] |>
        solve() |>
        diag() |>
        sqrt()) |> round(4)
