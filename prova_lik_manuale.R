rm(list = ls())
devtools::load_all()

# esempio 1: tutti design bivariati, solo baseline
des <- ncrr.design(smoke.alarm)
des <- subset(des,
              which(sapply(des$design, \(d) 0 %in% d)))

manlik3b <- function(params, theta, gamma, K, design) {
  # param in ordine: alpha, beta, mu0, rho, sigma20, sigma2
  alpha <- params[1:K]
  beta <- params[(K + 1):(2 * K)]
  mu0 <- params[2 * K + 1]
  rho <- 0.5#params[2 * K + 2]
  sigma20 <- exp(params[2 * K + 3])
  sigma2 <- exp(params[(2 * K + 4):length(params)])
  lik.vec <- sapply(seq_along(design), \(i) {
    d <- design[[i]][-1]
    mu <- c(mu0, alpha[d] + beta[d] * mu0)
    # funzione per selezionare la varianza giusta nel caso in cui sigma fosse
    # più corta del numero di parametri (ad es vcov.type = "achana")
    ds <- \(d) 1 + (d - 1) %% length(sigma2)
    if (length(d) == 1) {
      sigma <- matrix(c(sigma20 + gamma[[i]][1],
                        beta[d] * sigma20,
                        beta[d] * sigma20,
                        beta[d]^2 * sigma20 + sigma2[ds(d)] +
                          gamma[[i]][4]),
                      nrow = 2, ncol = 2)
    } else {
      sigma <- sigma20 * tcrossprod(c(1, beta[d]))
      sigma[2, 3] <- sigma[3, 2] <- sigma[3, 2] + rho * sqrt(prod(sigma2[ds(d)]))
      sigma[1, 1] <- sigma[1, 1] + gamma[[i]][1, 1]
      sigma[2, 2] <- sigma[2, 2] + sigma2[ds(d[1])] + gamma[[i]][2, 2]
      sigma[3, 3] <- sigma[3, 3] + sigma2[ds(d[2])] + gamma[[i]][3, 3]
      print(sigma)
    }
    sigmainv <- try(solve(sigma))
    if (inherits(sigmainv, "try-error")) {
      #browser()
      return(-Inf)
    }
    ss <- (\(x) as.vector(-0.5 * (x - mu) %*% sigmainv %*% (x - mu)))(theta[[i]])
    -0.5 * as.vector(determinant(sigma)$modulus) + ss
  })
  print(lik.vec)
  sum(lik.vec)
}

K <- 5
init <- getInitial(des, vcov.type = "achana", seed = 32)
# aggiungo il parametro rho per compatibilità con la funzione appena creata
init2b <- array(NA, length(init) + 1)
i <- j <- 1
while (i <= length(init)) {
  init2b[j] <- init[i]
  if (j == 2 * K + 2) {
    init2b[j] <- 0.5
  } else {
    i <- i + 1
  }
  j <- j + 1
}

llik <- get.llik.from.design(des, vcov.type = "achana", echo = 0)

theta <- crr.get.theta(des, raw = TRUE)
gamma <- crr.get.Gamma(des, raw = TRUE)

llik(init, theta, gamma)
manlik3b(init2b, theta, gamma, K, des$design)

nlminb(init2b, manlik3b, theta = theta, gamma = gamma, K = K, design = des$design)
nlminb(init, llik, y = theta, Gamma = gamma)

# a meno delle costanti moltiplicative le due funzioni mi sembrano equivalenti

# esempio 2: effetti normali, gruppi da 3, no baseline

des2 <- ncrr.design(morphine)
des2 <- subset(des2, which(sapply(des2$design, \(d) 0 %in% d)))

K <- 3
llik <- get.llik.from.design(des2, vcov.type = "achana", echo = 0)
init <- getInitial(des2, vcov.type = "achana", seed = 32)
# stavolta parto dalla massima verosimiglianza (si può anche togliere questa
# riga e proseguire con i punti di partenza casuali)
init <- optim(init, llik, method = "BFGS", control = list(fnscale = -1))$par
# aggiungo il parametro rho per compatibilità con la funzione appena creata
init2b <- array(NA, length(init) + 1)
i <- j <- 1
while (i <= length(init)) {
  init2b[j] <- init[i]
  if (j == 2 * K + 2) {
    init2b[j] <- 0.5
  } else {
    i <- i + 1
  }
  j <- j + 1
}

llik <- get.llik.from.design(des2, vcov.type = "achana", echo = 0)


theta <- crr.get.theta(des2, raw = TRUE)
gamma <- crr.get.Gamma(des2, raw = TRUE)

manlik3b(init2b, theta, gamma, K, des2$design)
llik(init, theta, gamma)

# il parametro rho alla posizione 2*K+2 non compare nella verosimiglianza e
# quindi va a caso
(pp2 <- nlminb(init, \(x) -llik(x, y = theta, Gamma = gamma)))

vcov.ncrr.design(des, pp2$par, llik = llik) |> diag() |> sqrt()

des2 <- subset(des2, sample(1:54, 100, TRUE))
des22 <- simulate.ncrr.design(des2, seed = 3)

# Metodo dei momenti

mom2b <- function(theta, gamma, K, design) {
  # INCOMPLETA
  strdes <- sapply(design, paste, collapse = "")
  des.idx <- Filter(\(x) length(x) > 1,
                    tapply(X = seq_along(des$design),
                           INDEX = sapply(design, paste, collapse = ""),
                           FUN = c))
  # come gestire design diversi in ogni caso?? Applicare direttamente i
  # principi della regressione con errori di misura non sembra fattibile. Che
  # non si debbano considerare i veri dati come un K-vettore per le medie e una
  # matrice di vcov K x K comuni a tutti gli studi?
  alpha <- beta <- sigma2 <- array(NA, K)
  mu0 <- rho <- sigma20 <- NA
  treat <- as.integer(substr(names(des.idx), 2, 2))
  for (i in seq_along(des.idx)) {
    th <- do.call(rbind, theta[des.idx[[i]]])
    mu <- colMeans(th)
    tau <- cov(th)
    alpha[t] <-  mu[2] - tau[1, 2] / tau[1, 1] * mu[1]
    beta[t] <-  tau[1, 2] / tau[1, 1]
    mu0 <-  mu[1] # ogni volta sovrascritto ?!?
    rho <-  NA
    sigma20 <-  tau[1, 1]
    sigma[t] <-  tau[2, 2] - tau[1, 2]^2 / tau[1, 1]
  }
}
