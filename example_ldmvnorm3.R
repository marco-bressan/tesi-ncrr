# trasformazione logaritmica dei parametri di varianza

rm(list = ls())
library(mvtnorm)
library(numDeriv)

.score1mat <- function(tt, mu, Sigma) {
  Scinv <- solve(Sigma)
  # score analitica per la media gaussiana V^(-1) sum(y_j - mu)
  scmu <- Scinv %*% (tt - mu)
  # score analitica per la varianza gaussiana
  scQ <- Scinv %*% tcrossprod(tt - mu) %*% Scinv
  scoreS <- -Scinv + 0.5 * diag(diag(Scinv)) + scQ - 0.5 * diag(diag(scQ))
  list(mu = scmu, Sigma = scoreS)
}

diagoffdiag <- function(x = 1, offx = 0, n = length(x)) {
  offx * (1 - diag(n)) + x * diag(n)
}

N <- 500L
J <- 5
set.seed(441)
# parametri reali
R <- diagoffdiag(1, .5, J)
diag(R) <- 1
Sigma <- diag(sqrt(1:J / 2)) %*% R %*% diag(sqrt(1:J / 2)) # varianze 0.5 1 1.5 2 2.5
mn <- c(2, 0, -1, 1, 2) # medie: mu0 = 1, beta = 0.0 -0.5  0.5  1.0

matrixcalc::is.positive.definite(R)

Y <- t(rmvnorm(N, mn, Sigma))
Yc <- Y - mn

veri <- c(mn[1], mn[-1] / mn[1], log(1:J / 2))


# riparametrizzo mu[1] = mu0, mu[i] = b[i-1] * mu0
# Sigma[i,j] = .5 * sqrt(s[i]) * sqrt(s[j]) se i != j altrimenti s[i]
# parm = [mu0   b[1] ... b[J-1]   s[1] ... s[J-1]]
ll5 <- function(parm) { # -loglik
  muc <- parm[1] * c(1, parm[2:J])
  S <- tcrossprod(sqrt(exp(parm[-(1:J)]))) * diagoffdiag(1, .5, J)
  sigmainv <- solve(S)
  ss <- apply(Y, 2, function(x) as.vector(-0.5 * (x - muc) %*% sigmainv %*% (x - muc)))
  sum(0.5*as.vector(determinant(S)$modulus)-ss)
}

sc5 <- function(parm) {
  R <- diagoffdiag(1, .5, J)
  muc <- parm[1] * c(1, parm[2:J])
  D <- diag(sqrt(exp(parm[-(1:J)])))
  S <- tcrossprod(diag(D)) * R
  sc <- apply(Y, 2, .score1mat, mu = muc, Sigma = S)
  dSig.ds <- mapply(\(s, j) {
    dDj <- D * 0
    dDj[j, j] <- .5 * 1/s
    (kronecker(D, dDj) + kronecker(dDj, D)) %*% c(R)
  }, diag(D), 1:J)
  ret <- vapply(sc, \(x) {
    c(mu0 = crossprod(x$mu, c(1, parm[2:J])),
      b = crossprod(x$mu, rbind(0, diag(parm[1], J-1))),
      s = crossprod(c(x$Sigma), dSig.ds))
  }, FUN.VALUE = numeric(2*J))
  ret <- rowSums(ret)
  ret[-(1:J)] <- ret[-(1:J)] * exp(parm[-(1:J)])
  -ret
}

sc52 <- function(parm) {
  R <- diagoffdiag(1, .5, J)
  muc <- parm[1] * c(1, parm[2:J])
  D <- diag(sqrt(exp(parm[-(1:J)])))
  S <- tcrossprod(diag(D)) * R
  sc <- apply(Y, 2, .score1mat, mu = muc, Sigma = S)
  dSig.ds <- sapply(1:J, \(j) {
    dSj <- D * 0
    for (t in 1:J)
      for (u in 1:J)
        dSj[t, u] <- ifelse(t == u && u == j,
                            1 / 2 / D[j, j],
                            ifelse(xor(u == j, t == j),
                                   S[t, u] / 2 / D[j, j]^2,
                                   0))
    c(dSj)
  })
  ret <- vapply(sc, \(x) {
    c(mu0 = crossprod(x$mu, c(1, parm[2:J])),
      b = crossprod(x$mu, rbind(0, diag(parm[1], J-1))),
      s = crossprod(c(x$Sigma), dSig.ds))
  }, FUN.VALUE = numeric(2*J))
  ret <- rowSums(ret)
  ret[-(1:J)] <- ret[-(1:J)] * exp(parm[-(1:J)])
  -ret
}

# implementazione di Sartori

sc5bis <- function(parm, rho = .5) {
  R <- diagoffdiag(1, rho, J)
  muc <- parm[1] * c(1, parm[2:J])
  Ycc <- Y - muc
  sigma <- exp(parm[-(1:J)])
  D <- diag(sqrt(sigma))
  S <- D %*% R %*% D
  invS <- solve(S)
  sc <- apply(Y, 2, .score1mat, mu = muc, Sigma = S)
  DinvSig.ds <- DSig.ds <- array(NA, c(J, J, J))
  # derivate delle entrate di Sigma
  for (t in 1:J)
    for (u in 1:J)
      for (j in 1:J)
        DSig.ds[t, u, j] <- ifelse(t == u && u == j,
                                    1 / 2 / D[j, j],
                                   ifelse(xor(u == j, t == j),
      # d/ds^2_j rho s_h s_j = rho s_h  * 1/(2 s_j) = S[h, j] / s_j * 1/(2 s_j)
                                          S[t, u] / 2 / sigma[j],
                                          0))
  # derivate delle entrate dell'inversa di Sigma
  for (r in 1:J)
    for (s in 1:J)
      for (j in 1:J)
        DinvSig.ds[r, s, j] <- invS[r, ] %*% DSig.ds[,, j] %*% invS[, s]
  # pezzi della derivata della verosimiglianza
  dlii.ds <- sapply(1:J, \(j) .5 * apply(Ycc, 2,
                                         \(y) t(y) %*% DinvSig.ds[,, j] %*% y))
  dli.ds <- sapply(1:J, \(j) .5 * sum(diag(invS %*% DSig.ds[,, j])))
  ret <- vapply(sc, \(x) {
    c(mu0 = crossprod(x$mu, c(1, parm[2:J])),
      b = crossprod(x$mu, rbind(0, diag(parm[1], J-1))))
  }, FUN.VALUE = numeric(J))
  ret <- c(rowSums(ret), rowSums(dli.ds + t(dlii.ds)))
  # derivata della trasformazione logaritmica
  ret[-(1:J)] <- ret[-(1:J)]
  -ret
}

# confronto con il gradiente numerico
cbind(sc5(veri), grad(ll5, veri))
cbind(sc52(veri), grad(ll5, veri))
cbind(sc5bis(veri), grad(ll5, veri))

# ottimizzazioni
(op5 <- nlminb(veri, ll5))
cbind(sc5bis(op5$par), grad(ll5, op5$par))
(op5 <- nlminb(veri, ll5, \(x) grad(ll5, x)))
(op52 <- nlminb(veri, ll5, \(x) sc5bis(x)))
(op53 <- nlminb(veri, ll5, \(x) sc5(x)))

sum((op5$par - veri)^2)
sum((op52$par - veri)^2)
sum((op53$par - veri)^2)

sc5bis(op5$par)
grad(ll5, op5$par)
