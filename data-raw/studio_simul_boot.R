#' ---
#' title: studio di simulazione
#' subtitle: dati di Achana
#' format:
#'   pdf:
#'     fig-width: 12
#'     fig-height: 10
#' ---
rm(list = ls())
#setwd("/home/marco/Documenti/tesi-ncrr/")
devtools::load_all()
library("likelihoodAsy")
#library("tesi.ncrr")

CONFLVL <- .95
NSIM <- 250
VCOVTYPE <- "achana"

#DIR <- "/home/marco/output-tesi"
#DIR <- "../.." # per il markdown
DIR <- "../output-boot/" # per l'esecuzione nel pacchetto
if (!dir.exists(DIR))
  dir.create(DIR)

#' # Simulazione basata sul problema di achana
des <- ncrr.design(smoke.alarm)
attr(des, "vcov.type") = "achana"

# stima MLE ottenuta sul dataset originale
simu.pars <- list(alpha = c(0.53118984013899, 1.0431973777787, 0.00434231242384523,
                            2.36407165618289, 2.66293182986318, 2.7339581049579),
                  beta = c(0.948918313002282, 1.02425107553256, 1.06604309779969,
                           0.240340161234799, 0.179383944934584, 0.179378353526596),
                  mu0 = 0.81098898311333,
                  sigma20 = 2.63212049308308,
                  sigma2 = 5.69469982077026) # stime MV dai dati originali
simu.des <- simulate(des, params = simu.pars,
                     vcov.type = VCOVTYPE, nsim = NSIM, seed = 212)
simu.pars.v <- tesi.ncrr:::crr.join.par(simu.pars) |> tesi.ncrr:::crr.transform.par()
gendat.fun <- function(data, theta) {
  theta <- crr.split.par(theta, length(data$treatments) - 1, transform = TRUE,
                         fixed = match.vcov.fixed(attr(data, "vcov.type")))
  # generazione dataset simulato (estraggo solo il design, non tutta la lista)
  suppressMessages(simulate.ncrr.design(data, params = theta, seed = 342)[[1]])
}

psi.fun <- function(theta, data) {
  theta[["beta5"]]
}



# bootstrap sul vero dataset
param <- match("beta5", names(simu.pars.v))
llik.fun <- get.llik.from.design2(des, vcov.type = VCOVTYPE, echo = 0,
                                  use.data = TRUE, stop.on.fail = FALSE)

#'
#' Si testa l'ipotesi che beta5 != 1. Da HMA:
#' eta_i = b + b_1 xi_i + e_i e_i ~ N (0, s^2) , i = 1,...,N
#' Usually, the inferential interest is in the
#' parameter β_1 associated with the underlying risk. If β_1 = 1, then if the
#' control risk increases by a certain amount, the treatment group risk
#' increases by the same amount. Thus, interesting cases are usually those
#' where β1 deviates from 1.
#'
boot.rp <- crr.boot(des, rp.stat, R = 500,
                    ran.gen = gendat.fun, mle = simu.pars.v, parallel = TRUE,
                    seed = c(1998135100L, 2044097286L, 1091132551L,
                             966088075L, 1553350452L, 1303502678L),
                    psi0 = 1, init = simu.pars.v, param = param, exact = TRUE,
                    retain.data = TRUE)
# controllo correttezza risultati
args <- c(list(des), attributes(boot.rp$t0), boot.rp$.dots)
debugonce(rp.stat)
do.call(rp.stat, args)
debugonce(rp.stat)
for (i in sample.int(1000, 10)) {
  args <- c(lapply(attributes(boot.rp)[c( "data", "theta.hat", "J" )], \(x) x[[i]]),
            boot.rp$.dots)
  names(args)[1] <- "dati.gen"
  args[["psi0"]] <- 2
  do.call(rp.stat, args)
}
# intervalli
bgrid <- seq(-10, 20, length.out = 50)
boot.stat <- Filter(is.finite, boot.rp$t)
density(boot.stat) |> plot()
abline(v = boot.rp$t0)
r.val <- lapply(bgrid, \(b) rp.stat(des, psi0 = b, init = simu.pars.v,
                                    param = param, exact = FALSE))
sigb2.val <- sapply(r.val, function(x) mean(boot.stat <= x))
sigb2.val <- clamp(sigb2.val, eps = 1e-8)
sm1 <- smooth.spline(qnorm(sigb2.val), bgrid)
# intervalli di confidenza
predict(sm1, qnorm(c(.975, .5, .025)))[["y"]]
#| eval: true
#| fig-cap: "Intervalli di confidenza bootstrap"
plot(sm1$y, sm1$x, type = "l", main = "Radice con segno del log-RV",
     xlab = expression(beta[paste(5, ",", i)]),
     ylab = expression(r[P](beta[paste(5, ",", i)])))
points(bgrid, qnorm(sigb2.val), pch = 20)
abline(h = qnorm(c(0.025, 0.975)), lty = 3)
abline(v = 1, lty = 2)

#versione automatica
boot.rp.ci1 <- crr.boot.ci(boot.rp, psi.grid = bgrid, within = FALSE, exact = FALSE, statistic = rp.stat)
boot.rp.ci2 <- crr.boot.ci(boot.rp, psi.grid = bgrid, within = NA, statistic = rp.stat, exact = TRUE,
                           parallel = TRUE)
# within = NA, exact = T :  700 s
# within = NA, exact = NA:  862 s
# within = T , exact = F : 1659 s

boot.rp.ci <- crr.boot.ci(boot.rp, psi.grid = bgrid[seq_along(bgrid)%%5 == 0],
                          statistic = rp.stat,
                          within = TRUE, exact = FALSE, parallel = TRUE)
#| eval: false
if (FALSE) {
  # differenza di performance tra la versione approssimata e quella esatta
  # della statistica rp (ci mette un po')
  microbenchmark::microbenchmark(
    exact = sapply(samp, \(x) rp.stat(des, x, init = simu.pars.v, param = param)),
    not_exact = sapply(samp, \(x) rp.stat(des, x, init = simu.pars.v, param = param,
                                          exact = FALSE)),
    times = 5,
    setup = {
      samp <- sample(bgrid, 10)
    }
  )
  #| eval: false
  # simulazione
  lik.vals <- lik.vals2 <- numeric(NSIM)
  init <- getInitial(simu.des[[1]], vcov.type = VCOVTYPE)
  par.J <- par.J2 <- array(NA, c(length(init), length(init), NSIM))
  par.h0 <- par.stime <- par.stime2 <- par.sd <- par.sd2 <- matrix(NA, length(init), NSIM)
  psi.rs <- psi.stime <- psi.sd <- numeric(NSIM)


  for (k in 1:NSIM) {
    message(sprintf("%.2f%%\r", k / NSIM * 100))
    llik.fun <- get.llik.from.design(simu.des[[k]], vcov.type = VCOVTYPE, echo = 0,
                                     use.data = TRUE)
    # ------ OTTIMIZZAZIONE ----
    opt1 <- crr.rstar(simu.des[[k]], thetainit = init, floglik = llik.fun,
                      fpsi = psi.fun,  psival = psi.fun(init),
                      datagen = gendat.fun, ronly = TRUE)
    unclass(opt1)
    boot.rp <- crr.boot(simu.des[[k]], rp.stat, R = 30, sim = "parametric",
                        ran.gen = gendat.fun, mle = opt1$theta.hat, parallel = TRUE,
                        seed = c(1998135100L, 2044097286L, 1091132551L, 966088075L, 1553350452L,
                                 1303502678L),
                        psi0 = simu.pars.v[param], init = opt1$theta.hat, param = param)

    saveRDS(boot.rp, file = file.path(DIR, paste0("boot_", k, "_", as.integer(Sys.time()), ".rds")))
    # ------ REGISTRAZIONE RISULTATI ----
    ## lik.vals[k] <- opt1$value
    ## par.stime[, k] <- tesi.ncrr:::crr.transform.par(opt1$par, inverse = TRUE)
    ## par.J[,,k] <- opt1$hessian
    ## par.sd[, k] <- sqrt(diag(solve(opt1$hessian)))
    ## if (inherits(opt2, "try-error")) next
    ## lik.vals2[k] <- llik.fun(opt2$theta.hat, simu.des[[k]])
    ## par.stime2[, k] <- tesi.ncrr:::crr.transform.par(opt2$theta.hat, inverse = TRUE)
    ## par.J2[,,k] <- opt2$info.hat
    ## par.sd2[, k] <- opt2$se.theta.hat
    ## par.h0[, k] <- opt2$theta.hyp
    ## # stime di psi
    ## psi.rs[k] <- opt2$rs
    ## psi.stime[k] <- opt2$psi.hat
    ## psi.sd[k] <- opt2$se.psi.hat
  }

  dimnames(par.h0) <- dimnames(par.stime) <- dimnames(par.stime2) <-
    dimnames(par.sd) <- dimnames(par.sd2) <-
    list(pars = names(init), repl = seq_along(simu.des))
  # fine simulazione
  save.image(file.path(DIR, "sim1provv"))

}
