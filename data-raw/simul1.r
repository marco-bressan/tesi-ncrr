#' ---
#' title: "prova dell'ottimizzatore per la network crr"
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
#' I dati sono gli stessi del primo esempio in Noma et al. sulla
#' cessazione dell'abitudine al fumo tramite 3 diverse tipologie di
#' assistenza.
rm(list = ls())
#| output: false
devtools::load_all(".")

smoking$tik <- with(smoking, log(rik) - log(nik - rik))
names(smoking)

# specifico il design della meta-analisi
des <- ncrr.design(smoking)
str(des)
des

#' # Riproduzione degli esempi giocattolo contenuti nel documento `verosim1.pdf`
#'
#' ## Due studi con design {0; 1} e {0; 2}
#'

# definisco il design `toy1` prendendo gli studi n. 4 e 6, che hanno design
# diversi, compatibili con l'esempio giocattolo
toy1 <- subset(des, c(4, 6))
init1 <- getInitial(toy1, transform = FALSE)
# creo la funzione di ottimizzazione a partire dal design specificato
fn1 <- get.llik.from.design(toy1, echo = 3, transform = FALSE)
#| output: false
# ottimizzazione vincolata
mv1 <- optim(init1, \(x) -fn1(x),
             lower = attr(init1, "lower"),
             upper = attr(init1, "upper"),
             method = "L-BFGS-B",
             control = list(fnscale = 1e-10, factr = 1, maxit = 1e6))

#'
#'
#' $\sigma^2_{01}$ e $\sigma^2_{02}$ invece risultano negativi!!
#'

# controllo l'esito dell'ottimizzazione
source("data-raw/diff_alpha1.R")
with(mv1, cat("Esito: ", convergence, " - ",
              if (is.null(message)) "OK" else message, "\n"))


# confronto i parametri calcolati con l'ottimizzatore...
(param1 <- crr.split.par(mv1$par, 2))

# ...con quelli in forma chiusa
(alphahat <- do.call(alpha.cf1, append(param1[c("sigma2", "beta", "sigma20", "mu0")],
                                      list(design = toy1)))
)# alpha calcolato in forma chiusa
(sigmahat <- do.call(sigma.cf1, append(param1[c("beta", "sigma20", "mu0")],
                                      list(design = toy1)))
) # sigma2 calcolato in forma chiusa
#'
#' Effettivamente, però, la forma della derivata (che ho calcolato con Maxima)
#' sembrerebbe essere proprio così
#'
#' \begin{align*}
#' \sigma_{01} & =-\frac{(\gamma_{11}+\beta_{01}^{2}\gamma_{10})\sigma_{0}+
#'   \gamma_{10}\gamma_{11}}{\sigma_{0}+\gamma_{10}}\\
#' \sigma_{02} & =-\frac{(\gamma_{22}+\beta_{02}^{2}\gamma_{20})\sigma_{0}+
#'   \gamma_{20}\gamma_{22}}{\sigma_{0}+\gamma_{20}}
#' \end{align*}
#'
#'
#' ::: {.callout-tip}
#'
#' ## **C'è un'apparente contraddizione!**
#'
#' Teoricamente $\sigma_{01}$ dovrebbe essere un parametro di varianza, quindi
#' positivo; ma nell'economia del modello, ossia se si considera semplicemente
#' la distribuzione di $\boldsymbol{\hat\theta}$, è sufficiente che
#' $\beta_{01}^{2}\sigma_{0}^{2}+\sigma_{01}^{2} > 0$ perchè la matrice di
#' Var-cov della normale sia sensata.
#'
#' Per cui, il modello "vuole andare" verso $\sigma$ negativi perchè tanto la
#' densità di $\boldsymbol{\hat\theta}$ è comunque ben definita, mentre invece
#' noi sappiamo che non lo è per via del modello statistico che sottende
#' questa distribuzione.
#'
#' Questo spiegherebbe anche il motivo per cui l'ottimizzatore s'incastra su 0.
#'
#' **Come si può agire in questo caso?**
#' :::
#'
#' Provo a vedere se le stime cambiano fissando i parametri
#'

# fisso parametro alpha
fixpar <- list(alpha = alphahat)
# disabilito stampa dei passaggi intermedi con `echo = 0`
fn1 <- get.llik.from.design(toy1, echo = 0, transform = FALSE)
# occhio che ogni volta bisogna fare in modo di togliere i parametri fissati!
mv1f <- optim(crr.remove.par(init1, names(fixpar)), \(x) -fn1(x, fixed = fixpar),
              lower = crr.remove.par(attr(init1, "lower"), names(fixpar)),
              upper = crr.remove.par(attr(init1, "upper"), names(fixpar)),
              method = "L-BFGS-B",
              control = list(fnscale = 1e-10, factr = 1, maxit = 1e6))

# confronto tra i parametri da stima non vincolata e vincolata rispettivamente
all.equal(crr.split.par(mv1$par, 2),
          crr.split.par(c(fixpar[[1]], mv1f$par), 2))

#' La differenza è comunque minima.

c(valore_ottimo = (mv1$value),
  valore_ottimo_vincolato = mv1f$value,
  differenza_pct = 100 * (mv1$value - mv1f$value)/mv1$value)
#'
#' ### Proviamo con un problema semplificato (Achana et al.)

fn1 <- get.llik.from.design(toy1, echo = Inf, transform = FALSE, vcov.type = "achana")
init1a <- getInitial(toy1, vcov.type = "achana", transform = FALSE, seed = 1)
mv1f <- optim(init1a, \(x) -fn1(x),
              lower = attr(init1a, "lower"),
              upper = attr(init1a, "upper"),
              method = "L-BFGS-B",
              control = list(fnscale = 1e-10, factr = 1, maxit = 1e6))

mv1f

# inversa dell'hessiana (in II colonna i valori dei parametri)

optimHess(mv1f$par, \(x) -fn1(x)) |>
  solve() |>
  diag() |>
  sqrt() |>
  cbind(mv1f$par)

#'
#'
#'
#'
#' ## Due studi con design uguale {0; 1}
#'

# seleziono solo studi n. 4 e 5
toy2 <- subset(des, c(4, 5))
str(toy2)
init2 <- getInitial(toy2, transform = FALSE)
fn2 <- get.llik.from.design(toy2, echo = 1, transform = FALSE)
#| output: false
mv2 <- optim(init2, \(x) -fn2(x),
             lower = attr(init2, "lower"),
             upper = attr(init2, "upper"),
             method = "L-BFGS-B",
             control = list(fnscale = 1e-10, factr = 1, maxit = 1e6))

#| output: true
with(mv2, cat("Esito: ", convergence, " - ",
              if (is.null(message)) "OK" else message, "\n"))

# i parametri ottimizzati in un formato più leggibile...
param2 <- crr.split.par(mv2$par, 1)
param2

# ... da confrontare con le stime analitiche:

alphahat <- do.call(alpha.cf2, append(param2[c("sigma2", "beta", "sigma20", "mu0")],
                                      list(design = toy2)))
alphahat
sigmahat <- do.call(sigma.cf2, append(param2[c("alpha", "beta", "sigma20", "mu0")],
                                      list(design = toy2)))
sigmahat

## Al momento questo codice da' errore, devo capire perchè...
## #'
## #' Ottimizzazione vincolata
## #'
## fixpar <- list(alpha = alphahat)
## fn2 <- get.llik.from.design(toy2, echo = 0, transform = FALSE)
## mv2f <- optim(crr.remove.par(init2, "alpha"), \(x) -fn2(x, fixed = fixpar),
##               lower = attr(init2, "lower"),
##               upper = attr(init2, "upper"),
##               method = "L-BFGS-B",
##               control = list(fnscale = 1e-10, factr = 1, maxit = 1e6))

## # confronto i parametri
## all.equal(crr.split.par(mv2$par, 1),
##           c(fixpar, crr.split.par(mv2f$par, 1)))

## all.equal(crr.split.par(mv2$par, 1),
##           subst.params(crr.split.par(mv2f$par, 1), fixpar))

## # confronto il valore ottenuto della funzione obiettivo
## c(valore_ottimo = (mv2$value),
##   valore_ottimo_vincolato = mv2f$value,
##   differenza = mv2$value - mv2f$value)


#'
#' # Esperimenti più complessi
#'
#' **TODO**: Qui andranno implementate le derivate esplicite una volta che la struttura
#' generale del modello statistico sarà definita.
#'
#' ## Esperimento 1
#'
#' specifico il design e provo a
#' condurre un esperimento selezionando un gruppo di studi con design "semplice"
#' (2 trattamenti a confronto, di cui uno è il baseline)
#'

# prendo tutti gli studi con design {0; 1}
des1 <- subset(des, which(sapply(des$design, \(d) identical(c(0, 1), d))))
str(des1)
par.init <- getInitial(des1, transform = FALSE)

# calcolo esplicitamente i parametri della normale per hat(theta)
# prova per i mu0
do.call(crr.get.mu,
        append(list(des1),
               crr.split.par(param.mv0$par, np = 1)[c("alpha", "beta", "mu0")]))
# prova per i sigma
do.call(crr.get.sigma,
        append(list(des1),
               crr.split.par(param.mv0$par, np = 1)[c("beta", "sigma20", "rho", "sigma2")]))
do.call(crr.get.sigma,
        append(list(des1, type = "simplified"),
               crr.split.par(par.init, np = 1)[c("beta", "sigma20", "rho", "sigma2")]))
do.call(crr.get.sigma,
        append(list(des1, type = "achana"),
               crr.split.par(par.init, np = 1)[c("sigma20")]))

opt.fn <- get.llik.from.design(des1, echo = 3, transform = FALSE)
#| output: false
param.mv1 <- optim(par.init, \(x) -opt.fn(x),
                  lower = attr(par.init, "lower"),
                  upper = attr(par.init, "upper"),
                  method = "L-BFGS-B",
                  control = list(fnscale = 1e-10, factr = 1, maxit = 1e6))
#| output: true
# Stime MV in formato leggibile
crr.split.par(param.mv1$par, 1, FALSE)

#'
#' -------
#'
#'
#' Ottimizzazione sulle trasformate dei parametri (log per la varianza, Fisher
#' per la correlazione)
#'

#| output: false
opt.fn <- get.llik.from.design(des1, echo = 3, transform = TRUE)
param.mv12 <- optim(getInitial(des1, transform = TRUE), \(x) -opt.fn(x),
                    method = "Nelder-Mead",
                    control = list(fnscale = 1e-10, factr = 1, maxit = 1e6))
#| output: true
crr.split.par(param.mv12$par, 1, TRUE)
hess.mv12 <- optimHess(param.mv12$par, \(x) -opt.fn(x))
solve(hess.mv12[-5, ][, -5]) |> diag() |> sqrt()

# con matrice di achana
opt.fn <- get.llik.from.design(des1, echo = 3, transform = TRUE,
                               vcov.type = "achana")
param.mv12a <- optim(getInitial(des1, transform = TRUE, vcov.type = "achana"),
                    \(x) -opt.fn(x),
                    method = "Nelder-Mead",
                    control = list(fnscale = 1e-10, factr = 1, maxit = 1e6))
#| output: true
crr.split.par(param.mv12$par, 1, TRUE)
hess.mv12a <- optimHess(param.mv12a$par, \(x) -opt.fn(x))
solve(hess.mv12a) |> diag() |> sqrt()




#'
#' Ottimizzazione con random start
#'
# levo l'output
opt.fn <- get.llik.from.design(des1, echo = 0, transform = TRUE)
par2 <- getInitial(des1, seed = c(alpha = 3, beta = 4, sigma20 = 6, sigma2 = 9),
                   rep = 50, transform = TRUE)
for (i in 1:nrow(par2)) {
  mvcurr <- optim(par2[i, ], \(x) -opt.fn(x),
                  method = "Nelder-Mead",
                  control = list(fnscale = 1e-10, factr = 1, maxit = 1e6))
  print(mvcurr$counts)
  par2[i, ] <- mvcurr$par
}

# stampa di statistiche descrittive per tutti i valori ottimizzati
crr.split.par(par2, transform = TRUE, np = 1)
crr.split.par(par2, transform = TRUE, np = 1) |>
  do.call(cbind, args = _) |>
  apply(2, \(x) c(min = min(x), max = max(x), mediana = median(x),
                  media = mean(x), dev.std = sd(x),
                  scarto.medio.assoluto.mediana = mean(abs(x - median(x))))) |>
  t() |>
  as.data.frame()

crr.split.par(par2, transform = TRUE, np = 1) |>
  do.call(cbind, args = _) |>
  apply(2, \(x) table(cut(x, breaks = 5)), simplify = FALSE)


#' ### Struttura di varianza-covarianza semplificata (sigma01=sigma02=...=sigma0)
#'
#'

#| output: false
opt.fn <- get.llik.from.design(des1, echo = 3, transform = TRUE,
                               vcov.type = "simplified")
param.mv12 <- optim(getInitial(des1, transform = TRUE), \(x) -opt.fn(x),
                    method = "Nelder-Mead",
                    control = list(fnscale = 1e-10, factr = 1, maxit = 1e6))
#| output: true
crr.split.par(param.mv12$par, 1, TRUE)
#'
#' Ottimizzazione con random start
#'
environment(opt.fn)$echo <- 0

par2 <- getInitial(des1, seed = c(alpha = 3, beta = 4, sigma20 = 6, sigma2 = 9),
                   rep = 50, transform = TRUE)
for (i in 1:nrow(par2)) {
  capture.output(type = "message",
                 mvcurr <- optim(par2[i, ], \(x) -opt.fn(x),
                                 method = "Nelder-Mead",
                                 control = list(fnscale = 1e-10, factr = 1, maxit = 1e6)))
  print(mvcurr$counts)
  par2[i, ] <- mvcurr$par
}

# stampa di statistiche descrittive per tutti i valori ottimizzati
crr.split.par(par2, transform = TRUE, np = 1) |>
  do.call(cbind, args = _) |>
  apply(2, \(x) c(min = min(x), max = max(x), mediana = median(x),
                  media = mean(x), dev.std = sd(x),
                  scarto.medio.assoluto.mediana = mean(abs(x - median(x))))) |>
  t() |>
  as.data.frame()

crr.split.par(par2, transform = TRUE, np = 1) |>
  do.call(cbind, args = _) |>
  apply(2, \(x) table(cut(x, breaks = 5)), simplify = FALSE)

#' ### Struttura di varianza-covarianza di Achana
#'
#'

#| output: false
opt.fn <- get.llik.from.design(des1, echo = 3, transform = TRUE,
                               vcov.type = "achana")
param.mv12 <- optim(getInitial(des1, transform = TRUE), \(x) -opt.fn(x),
                    method = "Nelder-Mead",
                    control = list(fnscale = 1e-10, factr = 1, maxit = 1e6))
#| output: true
crr.split.par(param.mv12$par, 1, TRUE)
#'
#' Ottimizzazione con random start
#'
environment(opt.fn)$echo <- 0

par2 <- getInitial(des1, seed = c(alpha = 3, beta = 4, sigma20 = 6, sigma2 = 9),
                   rep = 50, transform = TRUE)
for (i in 1:nrow(par2)) {
  capture.output(type = "message",
                 mvcurr <- optim(par2[i, ], \(x) -opt.fn(x),
                                 method = "Nelder-Mead",
                                 control = list(fnscale = 1e-10, factr = 1, maxit = 1e6)))
  print(mvcurr$counts)
  par2[i, ] <- mvcurr$par
}

# stampa di statistiche descrittive per tutti i valori ottimizzati
crr.split.par(par2, transform = TRUE, np = 1)[c("mu0", "beta", "alpha", "sigma20")] |>
  do.call(cbind, args = _) |>
  apply(2, \(x) c(min = min(x), max = max(x), mediana = median(x),
                  media = mean(x), dev.std = sd(x),
                  scarto.medio.assoluto.mediana = mean(abs(x - median(x))))) |>
  t() |>
  as.data.frame()

crr.split.par(par2, transform = TRUE, np = 1)[, c("mu0", "beta", "alpha", "sigma20")] |>
  do.call(cbind, args = _) |>
  apply(2, \(x) table(cut(x, breaks = 5)), simplify = FALSE)



#apply(par2, 1, \(x) crr.split.par(x, transform = TRUE, np = 1)$sigma20)
#apply(par2, 1, \(x) crr.split.par(x, transform = TRUE, np = 1)$rho)

#'
#' -----
#'
#' ## Esperimento 2
#'
#' Butto dentro tutti gli studi che hanno baseline 0, anche con design a tre trattamenti.
#'
#'

des2 <- subset(des, c(1, 3:4, 6:8))
str(des2)
par.init <- getInitial(des2, transform = TRUE)

# calcoli dei parametri della normale per hat(theta)
# prova per i mu0
do.call(crr.get.mu,
        append(list(des2),
               crr.split.par(par.init, np = 3, transform = TRUE)[c("alpha", "beta", "mu0")]))
# prova per i sigma
do.call(crr.get.sigma,
        append(list(des2),
               crr.split.par(par.init, np = 3, transform = TRUE)[c("beta", "sigma20", "rho", "sigma2")]))

#'
#' Conduco l'ottimizzazione sui parametri modificati (il log per i parametri di varianza,
#' la trasformata di Fisher per la correlazione) tramite Nelder-Mead.
#'

opt.fn <- get.llik.from.design(des2, echo = 0, transform = TRUE)

#| output: false
param.mv2 <- optim(par.init, \(x) -opt.fn(x),
                   #lower = attr(par.init, "lower"),
                   #upper = attr(par.init, "upper"),
                   method = "Nelder-Mead",
                   control = list(fnscale = 1e-10, factr = 1, maxit = 1e6))
#| output: true
# stime MV in formato leggibile
crr.split.par(param.mv2$par, 3, TRUE)

#'
#'
#'
#' Provo anche con il Simulated Annealing, come alternativa lenta ma che dovrebbe
#' fornire un risultato più robusto, non avendo a disposizione le derivate esplicite.
#'
#'

#| output: false
#| eval: false
sann.mv2 <- optim(par.init, \(x) -opt.fn(x),
                  #lower = attr(par.init, "lower"),
                  #upper = attr(par.init, "upper"),
                  method = "SANN",
                  control = list(fnscale = 1e-10, factr = 1, maxit = 1e6,
                                 temp = 10, tmax = 10, trace = 1))
save(sann.mv2, file = "sann_optim.rda")

#' Diagnostiche di confronto tra Nelder-mead e Simulated Annealing
#'

load("sann_optim.rda")

# confronto dei parametri
all.equal(crr.split.par(param.mv2$par, 3, TRUE),
          crr.split.par(sann.mv2$par, 3, TRUE))

# confronto dei valori della f.o.
c(valore_ottimo_NELMEAD = (param.mv2$value),
  valore_ottimo_SANN = sann.mv2$value,
  differenza = param.mv2$value - sann.mv2$value)

#'
#' Facendolo invece partire dall'ottimo individuato da Nelder-Mead:
#'
#| eval: false
sann.mv22 <- optim(param.mv2$par, \(x) -opt.fn(x),
                   #lower = attr(par.init, "lower"),
                   #upper = attr(par.init, "upper"),
                   method = "SANN",
                   control = list(fnscale = 1e-10, factr = 1, maxit = 1e6,
                                  temp = 10, tmax = 10, trace = 1))
save(sann.mv22, file = "sann_optim2.rda")

#' Confronti

load("../sann_optim2.rda")

# confronto dei parametri
all.equal(crr.split.par(param.mv2$par, 3, TRUE),
          crr.split.par(sann.mv22$par, 3, TRUE))

# confronto valori della f.o.
c(valore_ottimo_NELMEAD = (param.mv2$value),
  valore_ottimo_SANN = sann.mv22$value,
  differenza = param.mv2$value - sann.mv22$value)

#'
#' Nelder-Mead con punti di partenza casuali
#'

#| output: false
par22 <- getInitial(des2, seed = c(alpha = 3, beta = 4, sigma20 = 6, sigma2 = 9),
                   rep = 50, transform = TRUE)
for (i in 1:nrow(par22)) {
  mvcurr <- optim(par22[i, ], \(x) -opt.fn(x),
                  method = "Nelder-Mead",
                  control = list(fnscale = 1e-10, factr = 1, maxit = 1e6))
  message("Random start: ", i, "; n. iter: ", mvcurr$counts[1])
  par22[i, ] <- mvcurr$par
}

#| output: true
# per ogni random start calcolo alcune statistiche delle stime MV ottenute
# (sui parametri ritrasformati in scala originale)
crr.split.par(par22, 3, transform = TRUE) |>
  do.call(cbind, args = _) |>
  apply(2, \(x) c(min = min(x), max = max(x), mediana = median(x),
                  media = mean(x), dev.std = sd(x),
                  scarto.medio.assoluto.mediana = mean(abs(x - median(x))))) |>
  t() |>
  as.data.frame()

#' non un grande risultato...
