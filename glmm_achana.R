rm(list = ls())
library(lme4)

load("data/morphine.rda")

# versione di achana??

# In questa parte del codice provo a ricreare un modello lineare analogo a
# quello presentato nell'eq 6 del paper di achana, prendendo a modello la
# specificazione fatta nelle eqq 4-5.

trtmat <- model.matrix(~ 0 + treatment, morphine)
colnames(trtmat) <- levels(morphine$treatment)

# medie generali dei trattamenti baseline. Il primo elemento del vettore
# corrisponde alla media generale \overline\mu = 45.26 menzionata nel paper a
# pag. 759
mumb <- with(morphine[morphine$is.baseline == 1, ], tapply(mik, treatment, mean))

# ho costruito la matrice del disegno su modello della X contenuta nell'eq. 5;
# tuttavia, il problema di questa specificazione risiede nel fatto che tutte le
# \mu sono quantità aleatorie, essendo in un contesto di modello lineare con
# errore di misura. Dunque ho definito un modello lineare in cui la variabile
# esplicativa è il rischio osservato all'interno del baseline specifico per
# studio; questo è equivalente al modello di Achana fintanto che il baseline
# specifico corrisponde effettivamente al placebo, mentre per gli studi che non
# includono il placebo ho usato il rischio osservato per il baseline dello
# studio.
X <- tapply(morphine[, c("treatment", "mik")], morphine$study.id, \(dd) {
  n <- nrow(dd) - 1
  t <- as.integer(dd$treatment)
  t <- cbind(t[-1], t[1]) - 1
  r <- matrix(0, n, length(levels(morphine$treatment)) - 1)
  r[(1L:n) + n * (t[, 1] - 1L)] <- 1 # nnet::class.ind
  r[(t[, 2] != 0) * ((1L:n) + n * (t[, 2] - 1L))] <- -1
  rbind(0, cbind(r, r * (dd$mik[1] - mumb[1])))
}) |>
  do.call(what = rbind) |>
  as.data.frame()
muib <- with(morphine, tapply(mik, study.id, \(x) rep(x[1], length(x)))) |>
  unlist()
# i nomi corrispondono con la notazione usata nel paper, ma gli indici sono
# invertiti (3=paracetamol, 2 = nsaid, 1 = cox2)
names(X) <- c("d01", "d02", "d03", "beta01", "beta02", "beta03")

# Nel modello lineare tolgo le intercette e applico un offset pari al valore
# osservato del rischio baseline specifico per lo studio, in modo che la
# specificazione del modello risulti analoga a quella dell'eq. 6.
summary(fe <- lm(morphine$mik ~ 0 + offset(muib) + ., data = X))
# In questa seconda specificazione includo anche un peso inversamente
# proporzionale alla varianza within-study osservata, su modello della classica
# Weighted Least Squares menzionata nell'Handbook of meta-analysis e nei vari
# paper di Guolo.
summary(few <- lm(morphine$mik ~ 0 + offset(muib) + ., data = X,
                  weights = sqrt(morphine$nik) / morphine$sik))
# confrontando i beta del modello pesato con il C1 del paper (tab III) si nota che i
# risultati non sono troppo diversi, ma gli intervalli sono più stretti, il che
# è comprensibile dato che l'errore di misura non viene modellato in alcun modo
# e i muib vengono assunti osservati.
cbind(est = coef(few), confint(few))

# ho provato anche un modello ad effetti casuali ma non sono sicuro di come
# specificarlo correttamente... Infatti questo stima molti beta ~= 0

# stima simil-Salanti et al.
summary(re <- lmer(morphine$mik ~ 0 + offset(muib) + d01 + d02 + d03 +
                     (beta01 | morphine$study.id) + (beta02 | morphine$study.id) +
                     (beta03 | morphine$study.id),
                   data = X))
# stima solo intercette casuali (non converge nemmeno...)
summary(re2 <- lmer(morphine$mik ~ 0 + offset(muib) + . + (1 | morphine$study.id),
                    data = X))



# implementazione in jags

# per comodità, ho costruito il dataset copiando "a mano" i valori riportati
# nella tabella in appendice

jags_data_raw <- read.csv("data-raw/csv/jags data.csv", colClasses = "numeric")

jags_data <- with(
  jags_data_raw,
  list(
    ns = nrow(jags_data_raw),
    nt = 4,
    id = jags_data_raw$id,
    na = jags_data_raw$na,
    t = as.matrix(jags_data_raw[, grep("^t", names(jags_data_raw))]),
    y = as.matrix(jags_data_raw[, grep("^y", names(jags_data_raw))]),
    se = as.matrix(jags_data_raw[, grep("^se", names(jags_data_raw))]),
    bs_mean = 45.26
  ))

# inizializzazione
inits <- list(
  d = c(NA, 0, 0, 0), sd = 1,
  mu = rep(0, nrow(jags_data_raw)),
  delta = matrix(as.integer(is.na(jags_data$t) | jags_data$t == 1), ncol = 3),
  y = matrix(as.integer(is.na(jags_data$y) &
                          !is.na(jags_data$t) & jags_data$t == 1),
             ncol = 3),
  sdmu = 1, mu_mean = 0, beta = c(NA, 0, 0, 0)
)
inits$delta[inits$delta == 1] <- NA
inits$y[inits$y == 0] <- NA

library(rjags)

# fitting del modello JAGS (dovrebbe essere equivalente a quello riportato nel
# paper, l'ho fatto revisionare da chatGPT per la conversione da WinBUGS a
# JAGS)
fit <- jags.model(file = "achana-model-2.jags",
                  data = jags_data,
                  inits = inits,
                  n.chains = 3, n.adapt = 0)
adapt(fit, 1000) # adattamento OK
update(fit, n.iter = 30000) # burnin

beta.samp <- coda.samples(fit, "beta", n.iter = 70000)
summary(beta.samp)
plot(beta.samp)
coda::effectiveSize(beta.samp)

# come si vede dal summary i parametri beta 2 e 3 sono diversi da quelli
# riportati nei risultati del paper
