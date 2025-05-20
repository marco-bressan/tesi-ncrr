#' ---
#' title: "prova dell'ottimizzatore per la network crr"
#' format:
#'   pdf:
#'     fig-width: 8
#'     fig-height: 6
#' execute:
#'   eval: true
#' ---
#'
#' # Importazione dei dati
#'
#' I dati sono gli stessi del primo esempio in Noma et al. sulla
#' cessazione dell'abitudine al fumo tramite 3 diverse tipologie di
#' assistenza.

#| output: false
devtools::load_all(".")

data("smoking", package = "tesi.ncrr")
ls()
smoking$tik <- with(smoking, log(rik) - log(nik - rik))
names(smoking)

des <- ncrr.design(smoking)
str(des)
des

#' # Prima prova
#'
#' codice "legacy" in cui provo un singolo studio (il 3 : A (baseline) vs. B)

theta.oss <- smoking[smoking$study.id == 3, "tik"]
Gamma <- crr.vcov.within(smoking[smoking$study.id == 3, "rik"],
                         smoking[smoking$study.id == 3, "nik"])
#| output: false
param.mv0 <- optim(c(0, 1, 0, 1, 0, 1), \(x) -llik1(x, theta.oss, Gamma),
                  lower = c(-Inf, -Inf, -Inf, 1e-10, -1, 1e-10),
                  upper = c(Inf, Inf, Inf, Inf, 1, Inf),
                  method = "L-BFGS-B",
                  control = list(fnscale = 1e-10, factr = 1, maxit = 1e6))


#' ottengo stime dei parametri sensate solo per i parametri di posizione, mentre
#' le varianze e la correlazione si "incastrano su 0"

param.mv0


#' # Riproduzione di esempi giocattolo contenuti nel documento `verosim1`
#'
#' ## Due studi con design {0; 1} e {0; 2}
#'

toy1 <- subset(des, c(4, 6))
str(toy1)
init1 <- getInitial(toy1, transform = FALSE)
fn1 <- get.llik.from.design(toy1, echo = 0, transform = FALSE)
#| output: false
mv1 <- optim(init1, \(x) -fn1(x),
             lower = attr(init1, "lower"),
             upper = attr(init1, "upper"),
             method = "L-BFGS-B",
             control = list(fnscale = 1e-10, factr = 1, maxit = 1e6))


#'
#'
#' $\sigma2_{01}$ e $\sigma2_{02}$ invece risultano negativi!!
#'

source("diff_alpha1.R")
with(mv1, cat("Esito: ", convergence, " - ",
              if (is.null(message)) "OK" else message, "\n"))


param1 <- crr.split.par(mv1$par, 2)
param1
alphahat <- do.call(alpha.cf1, append(param1[c("sigma2", "beta", "sigma20", "mu0")],
                                      list(design = toy1)))
alphahat
sigmahat <- do.call(sigma.cf1, append(param1[c("beta", "sigma20", "mu0")],
                                      list(design = toy1)))
sigmahat
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
#' ::: {.callout-warninig}
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
#' Si può provare a vedere se le stime cambiano fissando i parametri
#'

fixpar <- list(alpha = alphahat)
fn1 <- get.llik.from.design(toy1, echo = 0, transform = FALSE)
mv1f <- optim(init1, \(x) -fn1(x, fixed = fixpar),
              lower = attr(init1, "lower"),
              upper = attr(init1, "upper"),
              method = "L-BFGS-B",
              control = list(fnscale = 1e-10, factr = 1, maxit = 1e6))


all.equal(crr.split.par(mv1$par, 2),
          subst.params(crr.split.par(mv1f$par, 2), fixpar))

c(valore_ottimo = (mv1$value),
  valore_ottimo_vincolato = mv1f$value,
  differenza = mv1$value - mv1f$value)

#'
#'
#' ## Due studi con design uguale {0; 1}
#'

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

with(mv2, cat("Esito: ", convergence, " - ",
              if (is.null(message)) "OK" else message, "\n"))


param2 <- crr.split.par(mv2$par, 1)
param2
do.call(alpha.cf2, append(param2[c("sigma2", "beta", "sigma20", "mu0")],
                          list(design = toy2)))
do.call(sigma.cf2, append(param2[c("alpha", "beta", "sigma20", "mu0")],
                          list(design = toy2)))

#'


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

des1 <- subset(des, which(sapply(des$design, \(d) identical(c(0, 1), d))))
str(des1)
par.init <- getInitial(des1, transform = FALSE)

# prova per i mu0
do.call(crr.get.mu,
        append(list(des1),
               crr.split.par(param.mv0$par, np = 1)[c("alpha", "beta", "mu0")]))
# prova per i sigma
do.call(crr.get.sigma,
        append(list(des1),
               crr.split.par(param.mv0$par, np = 1)[c("beta", "sigma20", "rho", "sigma2")]))
do.call(crr.get.sigma,
        append(list(des1),
               crr.split.par(par.init, np = 1)[c("beta", "sigma20", "rho", "sigma2")]))

opt.fn <- get.llik.from.design(des1, echo = 3, transform = FALSE)
#| output: false
param.mv1 <- optim(par.init, \(x) -opt.fn(x),
                  lower = attr(par.init, "lower"),
                  upper = attr(par.init, "upper"),
                  method = "L-BFGS-B",
                  control = list(fnscale = 1e-10, factr = 1, maxit = 1e6))

param.explicit <- crr.split.par(param.mv1$par, 1, FALSE)
# alfa esplicito: α_01=-(((t_10·β_01-t_11)·σ_0+(µ_0·β_01-t_11)·γ_10)/(σ_0+γ_10))

with(param.explicit,
     -((des1$theta[1] * beta - des1$theta[2]) * sigma20 +
       (mu0 * beta - des1$theta[2]) * des1$gamma[1]) / (sigma20 + des1$gamma[1]))


opt.fn <- get.llik.from.design(des1, echo = 3, transform = TRUE)
#| output: false
param.mv12 <- optim(getInitial(des1, transform = TRUE), \(x) -opt.fn(x),
                    method = "Nelder-Mead",
                    control = list(fnscale = 1e-10, factr = 1, maxit = 1e6))
#| output: true
param.mv12
#'
#' con random start
#'
environment(opt.fn)$echo <- 0

par2 <- getInitial(des1, seed = c(alpha = 3, beta = 4, sigma20 = 6, sigma2 = 9),
                   rep = 50, transform = TRUE)
for (i in 1:nrow(par2)) {
  mvcurr <- optim(par2[i, ], \(x) -opt.fn(x),
                  method = "Nelder-Mead",
                  control = list(fnscale = 1e-10, factr = 1, maxit = 1e6))
  print(mvcurr$counts)
  par2[i, ] <- mvcurr$par
}
# controllo
apply(par2, 1, \(x) crr.split.par(x, transform = TRUE, np = 1)$sigma20)
apply(par2, 1, \(x) crr.split.par(x, transform = TRUE, np = 1)$rho) # stime erratiche perchè non compare!


#'
#' -----
#'
#' ## Esperimento 2
#'
#' Butto dentro tutti gli studi che hanno baseline 0, anche con design a tre trattamenti.
#'
#'

des2 <- subset(des, c(1, 3:4, 6:8))
par.init <- getInitial(des2, transform = TRUE)
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
crr.split.par(param.mv2$par, 3, TRUE)
#'
#'
#' ::: {.callout-tip}
#'
#' ## **Integrazione delle derivate analitiche**
#'
#' Immagino che in questo caso la derivata debba tener conto della trasformazione
#' dei parametri, giusto?
#' :::
#'
#'
#'
#| output: false
par22 <- getInitial(des2, seed = c(alpha = 3, beta = 4, sigma20 = 6, sigma2 = 9),
                   rep = 50, transform = TRUE)
for (i in 1:nrow(par22)) {
  mvcurr <- optim(par22[i, ], \(x) -opt.fn(x),
                  method = "Nelder-Mead",
                  control = list(fnscale = 1e-10, factr = 1, maxit = 1e6))
  cat("Random start:", i, "; n. iter:", mvcurr$counts[1], "\n")
  par22[i, ] <- mvcurr$par
}

#| output: true
par22
