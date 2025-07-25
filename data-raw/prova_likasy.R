rm(list = ls())
devtools::load_all()
library(likelihoodAsy)


VCOVTYPE <- "achana"
#DIR <- "../.." # per il markdown
DIR <- ".." # per l'esecuzione nel pacchetto

#| warning: false
simu.pars <- list(alpha = c(0.53118984013899, 1.0431973777787, 0.00434231242384523,
                            2.36407165618289, 2.66293182986318, 2.7339581049579),
                  beta = c(0.948918313002282,
                           1.02425107553256, 1.06604309779969,
                           0.240340161234799, 0.179383944934584,
                           0.179378353526596),
                  mu0 = 0.81098898311333,
                  sigma20 = 2.63212049308308,
                  sigma2 = 5.69469982077026) # stime MV dai dati originali
simu.pars.v <- crr.join.par(simu.pars) |> crr.transform.par()

des <- ncrr.design(smoke.alarm, VCOVTYPE)
simu.des <- simulate(des, params = simu.pars, seed = 342)[[1]]
#data <- simulate(des, params = simu.pars)
llik <- get.llik.from.design(simu.des)
llik.fun <- function(theta, data) {
  y = crr.get.theta(data, raw = TRUE)
  Gamma = crr.get.Gamma(data, raw = TRUE)
  L <- llik(theta, y = y, Gamma = Gamma)
  # if (!is.finite(L)) {
  #   environment(llik)$echo <- 3
  #   browser()
  #   llik(theta, y = y, Gamma = Gamma)
  # }
  return(L)
}
gendat.fun <- function(theta, data) {
  theta <- crr.split.par(theta, length(data$treatments) - 1,
                         fixed = match.vcov.fixed(attr(data, "vcov.type")))
  suppressMessages(simulate(data, params = theta, seed = 342)[[1]]) # estraggo solo il design, non tutta la lista
}
psi.fun <- function(theta) {
  theta <- unlist(theta[grep("sigma", names(theta))])
  nm <- names(theta)
  theta <- .ptrans(theta, "sigma", inverse = TRUE)
  names(theta) <- nm
  L <- sum(theta)
  if (!is.finite(L)) browser()
  return(L)
}

options(warn = 2)
library(likelihoodAsy)
Rstar <- rstar(simu.des, thetainit = simu.pars.v, floglik = llik.fun,
               fpsi = psi.fun,  psival = psi.fun(simu.pars.v),
               datagen = gendat.fun, seed = 12, constr.opt = "alabama", R=250)
str(Rstar)
summary(Rstar)


# R = 250: seed = 544 ok, seed = 12 warning ma produce output
Rstar2 <- rstar(simu.des, thetainit = simu.pars.v, floglik = llik.fun,
                fpsi = psi.fun,  psival = psi.fun(simu.pars.v),
                datagen = gendat.fun, seed = 22, constr.opt = "solnp", R=250)
str(Rstar2)
summary(Rstar2)


save.image("prova_likasy1")

load("prova_likasy1")

# fare simulazione: uno dei sigma e qualche beta (beta5 che va peggio)
# confrontare stime di mv del pacchetto e mie
# verificare normalità di rstar

# rstarc
