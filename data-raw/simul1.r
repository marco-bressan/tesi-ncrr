# prova dell'ottimizzatore per la network crr

devtools::load_all(".")

data("smoking", package = "tesi.ncrr")
ls()
smoking$tik <- with(smoking, log(rik) - log(nik - rik))
names(smoking)

des <- ncrr.design(smoking)
str(des)
des

#' codice "legacy" in cui provo un singolo studio (il 3 : A (baseline) vs. B)
#' ottengo stime dei parametri sensate solo per i parametri di posizione, mentre
#' le varianze e la correlazione si "incastrano su 0"

theta.oss <- smoking[smoking$study.id == 3, "tik"]
Gamma <- crr.vcov.within(smoking[smoking$study.id == 3, "rik"],
                         smoking[smoking$study.id == 3, "nik"])
param.mv0 <- optim(c(0, 1, 0, 1, 0, 1), \(x) -llik1(x, theta.oss, Gamma),
                  lower = c(-Inf, -Inf, -Inf, 1e-10, -1, 1e-10),
                  upper = c(Inf, Inf, Inf, Inf, 1, Inf),
                  method = "L-BFGS-B",
                  control = list(fnscale = 1e-10, factr = 1, maxit = 1e6))

#'
#' # Esperimento 1
#'
#' specifico il design e provo a
#' condurre un esperimento selezionando un gruppo di studi con design "semplice"
#' (2 trattamenti a confronto, di cui uno è il baseline)
#'

des1 <- subset(des, which(sapply(des$design, \(d) identical(c(0, 1), d))))
str(des1)
par.init <- getInitial(des1)
# prova per i mu0
do.call(crr.get.mu,
        append(list(des1),
               crr.split.par(param.mv1$par, np = 1)[c("alpha", "beta", "mu0")]))
# prova per i sigma
do.call(crr.get.sigma,
        append(list(des1),
               crr.split.par(param.mv1$par, np = 1)[c("beta", "sigma20", "rho", "sigma2")]))
do.call(crr.get.sigma,
        append(list(des1),
               crr.split.par(par.init, np = 1)[c("beta", "sigma20", "rho", "sigma2")]))

opt.fn <- get.llik.from.design(des1)
param.mv1 <- optim(par.init, \(x) -opt.fn(x),
                  lower = attr(par.init, "lower"),
                  upper = attr(par.init, "upper"),
                  method = "L-BFGS-B",
                  control = list(fnscale = 1e-10, factr = 1, maxit = 1e6))



#'
#' -----
#'
#' # Esperimento 2
#'
#' il codice seguente non funziona ancora, devo sistemare il modo in cui
#' get.mu passa il vettore dei parametri a mean.baseline0
#'
#'

des2 <- subset(des, c(1, 3:4, 6:8))
par.init <- getInitial(des2)
# prova per i mu0
do.call(crr.get.mu,
        append(list(des2),
               crr.split.par(par.init, np = 3)[c("alpha", "beta", "mu0")]))
# prova per i sigma
do.call(crr.get.sigma,
        append(list(des2),
               crr.split.par(param.mv1$par, np = 1)[c("beta", "sigma20", "rho", "sigma2")]))
do.call(crr.get.sigma,
        append(list(des2),
               crr.split.par(par.init, np = 1)[c("beta", "sigma20", "rho", "sigma2")]))

opt.fn <- get.llik.from.design(des2)
param.mv1 <- optim(par.init, \(x) -opt.fn(x),
                  lower = attr(par.init, "lower"),
                  upper = attr(par.init, "upper"),
                  method = "L-BFGS-B",
                  control = list(fnscale = 1e-10, factr = 1, maxit = 1e6))
