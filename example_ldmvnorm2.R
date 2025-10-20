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

N <- 1000L
J <- 5
set.seed(2121)
# parametri reali
R <- sample(c(outer(c(-1, 0, 0, .35, 1), c(.25, .5, .75))), J^2, replace = TRUE) |>
  matrix(nrow = J)
R[upper.tri(R)] <- t(R)[upper.tri(R)]
diag(R) <- 1
Sigma <- diag(sqrt(1:J / 2)) %*% R %*% diag(sqrt(1:J / 2)) # varianze
C <- t(chol(Sigma))
mn <- 1:J # medie

matrixcalc::is.positive.definite(R)

# conversione in ltMatrices
prm <- C[lower.tri(C, diag = TRUE)]
lt <- ltMatrices(matrix(prm, ncol = 1L),
                   diag = TRUE, ### has diagonal elements
                   byrow = FALSE) ### prm is column-major
all.equal(C, as.array(lt)[,,1], check.attributes = FALSE, check.class = FALSE)
all.equal(Sigma, as.array(Tcrossprod(lt))[,,1], check.attributes = FALSE, check.class = FALSE)

# Vogliamo stimare C in max verosimiglianza
# We generate some data from N (0 , Σ) by first sampling from Z ~ N (0J , I ) and then
# computing Y = CZ + μ ∼ N_J (μ, CC^⊤)

Z <- matrix(rnorm(N * J), nrow = J)
Y <- Mult(lt, Z) + mn

## stime di massima verosimiglianza analitiche
(muhat <- rowMeans(Y))
(Shat <- var(t(Y)) * (N - 1) / N)

# stime numeriche per la mat. di vcov
Yc <- Y - rowMeans(Y) # centratura della Y

ll <- function(parm) {# -loglik
  muc <- parm[1:J]
  Cc <- try(chol(ss <- syMatrices(parm[-(1:J)], diag = TRUE, byrow = FALSE)),
            silent = TRUE)
  if (inherits(Cc, "try-error"))
    return(1e10)
  -ldmvnorm(obs = Y, mean = muc, chol = Cc)
}
sc <- function(parm) {# score
  muc <- parm[1:J]
  S <- as.array(syMatrices(parm[-(1:J)], diag = TRUE, byrow = FALSE))[, , 1]
  sc <- apply(Y, 2, .score1mat, mu = muc, Sigma = S)
  -rowSums(vapply(sc,
                  \(x) c(x[[1]], x[[2]][lower.tri(x[[2]], diag = TRUE)]),
                  numeric(length(parm))))
}

ML <- c(unname(muhat), Shat[lower.tri(Shat, diag = TRUE)])
ll(ML)
sc(ML)
start <- c(rep(0, J), diag(J)[lower.tri(diag(J), diag = TRUE)])
all.equal(numDeriv::grad(ll, start), sc(start),
          check.attributes = FALSE, check.class = FALSE)
llim <- c(rep(-Inf, J), -1 / diag(Inf, J)[lower.tri(diag(J), diag = TRUE)])
op <- optim(start, fn = ll, gr = sc, method = "L-BFGS-B",
             lower = llim, control = list(trace = FALSE))
all.equal(op$par, ML) ## OK
c(ll(op$par), ll(ML)) |> as.character()

## microbenchmark::microbenchmark(
##   score = optim(start, fn = ll, gr = sc, method = "L-BFGS-B",
##                 lower = llim, control = list(trace = FALSE)),
##   no_score = optim(start, fn = ll, method = "L-BFGS-B",
##                    lower = llim, control = list(trace = FALSE)),
##   times = 5
## )|>
##   (\(x) {print(x); x}) () |>
##                          plot()

# trasformazione logaritmica dei parametri di varianza

whichlog <- J * (J + 1) / 2 + 1 - cumsum(1:J)
ll4 <- function(parm) { # -loglik
  parm[J + whichlog] <- exp(parm[J + whichlog])
  muc <- parm[1:J]
  Cc <- try(chol(ss <- syMatrices(parm[-(1:J)], diag = TRUE, byrow = FALSE)),
            silent = TRUE)
  if (inherits(Cc, "try-error"))
    return(1e10)
  -ldmvnorm(obs = Y, mean = muc, chol = Cc)
}
sc4 <- function(parm) {
  parm[J + whichlog] <- exp(parm[J + whichlog])
  muc <- parm[1:J]
  S <- as.array(syMatrices(parm[-(1:J)], diag = TRUE, byrow = FALSE))[, , 1]
  sc <- apply(Y, 2, .score1mat, mu = muc, Sigma = S)
  ret <- rowSums(vapply(sc,
                        \(x) c(x[[1]], x[[2]][lower.tri(x[[2]], diag = TRUE)]),
                        numeric(length(parm))))
  ret[J + whichlog] <- ret[J + whichlog] * parm[J + whichlog]
  -ret
}
ML4 <- ML
ML4[J + whichlog] <- log(ML4[J + whichlog])
all.equal(ll4(ML4), ll(ML))
all.equal(sc4(ML4), sc(ML))
start4 <- start
start4[J + whichlog] <- log(start4[J + whichlog])
all.equal(ll4(start4), ll(start))
all.equal(sc4(start4), sc(start))

(op4 <- nlminb(start4, ll4, sc4))
(op42 <- nlminb(start4, ll4))
op4par <- op4$par
op4par[J + whichlog] <- exp(op4par[J + whichlog])
all.equal(op4par, op$par)

## microbenchmark::microbenchmark(
##   score = nlminb(start4, ll4, sc4),
##   no_score = nlminb(start4, ll4),
##   times = 5
## )|>
##   (\(x) {print(x); x}) () |>
##                          plot()


# trasformazione logaritmica dei parametri di varianza

# riparametrizzo mu[1] = mu0, mu[i] = b[i-1] * mu0
# Sigma[i,j] = .5 * sqrt(s[i]) * sqrt(s[j]) se i != j altrimenti s[i]
# parm = [mu0   b[1] ... b[J-1]   s[1] ... s[J-1]]
ll5 <- function(parm) { # -loglik
  muc <- parm[1] * c(1, parm[2:J])
  S <- tcrossprod(sqrt(exp(parm[-(1:J)]))) * diagoffdiag(1, .5, J)
  Cc <- try(chol(S), silent = TRUE)
  if (inherits(Cc, "try-error")) return(1e10)
  -ldmvnorm(obs = Y, mean = muc,
            chol = ltMatrices(t(Cc)[lower.tri(Cc, diag = TRUE)], diag = TRUE))
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

ML5 <- c(mu0 = ML[1], b = ML[2:J] / ML[1],
         s = optim(rep(1, J),
                   \(x) sum((Shat -
                               diagoffdiag(1, .5, J) * tcrossprod(sqrt(exp(x))))^2),
                   control = list(maxit = 2e6))$par)
MLbis <- c(ML[1:J], (diagoffdiag(1, .5, J) * tcrossprod(sqrt(exp(ML5[-(1:J)]))))[lower.tri(diag(J), diag = TRUE)])
all.equal(ll5(ML5), ll(MLbis)) # la funzione restituisce lo stesso valore?
start5 <- rnorm(length(ML5))

(op5 <- nlminb(start5, ll5, sc5))
(op52 <- nlminb(start5, ll5, control = list(eval.max = 30000, iter.max = 30000)))

cbind(sc5(start5), grad(ll5, start5))
cbind(sc5(op5$par), grad(ll5, op5$par))
cbind(sc5(op52$par), grad(ll5, op52$par))

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
  for (t in 1:J)
    for (u in 1:J)
      for (j in 1:J)
        DSig.ds[t, u, j] <- ifelse(t == u && u == j,
                                   .5 * 1 / D[j, j],
                                   ifelse(xor(u == j, t == j),
                                          S[t, u] / D[j, j],
                                          0))
  for (r in 1:J)
    for (s in 1:J)
      for (j in 1:J)
        DinvSig.ds[r, s, j] <- invS[r, ] %*% DSig.ds[,, j] %*% invS[, s]
  dlii.ds <- sapply(1:J, \(j) .5 * apply(Ycc, 2,
                                         \(y) t(y) %*% DinvSig.ds[,, j] %*% y))
  dli.ds <- sapply(1:J, \(j) .5 * sum(diag(invS %*% DSig.ds[,, j])))
  ret <- vapply(sc, \(x) {
    c(mu0 = crossprod(x$mu, c(1, parm[2:J])),
      b = crossprod(x$mu, rbind(0, diag(parm[1], J-1))))
  }, FUN.VALUE = numeric(J))
  ret <- c(rowSums(ret), rowSums(dli.ds + t(dlii.ds)))
  ret[-(1:J)] <- ret[-(1:J)] * exp(parm[-(1:J)])
  -ret
}

cbind(sc5bis(start5), grad(ll5, start5))

(op5 <- nlminb(start5, ll5))
(op5 <- nlminb(start5, ll5, \(x) sc5bis(x)))

sc5bis(op5$par)
