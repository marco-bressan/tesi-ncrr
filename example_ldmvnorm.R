### five observations
N <- 5000L
### dimension
J <- 4

devtools::load_all()
library(mvtnorm)
library(numDeriv)

### Cholesky factor
C <- ltMatrices(1/rnorm(J * (J + 1) / 2)^2, diag = TRUE)
### corresponding covariance matrix
S <- as.array(Tcrossprod(C))[,,1]

### lower and upper bounds, ie interval-censoring
obs <- mvtnorm::rmvnorm(N, sigma = S)

m <- matrix(rnorm(J * N), nrow = J)

### compare with pmvnorm
LL <- list(
  ldmvnorm = ldmvnorm(t(obs), chol = C, logLik = FALSE),
  lpmvnorm = apply(obs, 1, function(x) dmvnorm(x, sigma = S, log = TRUE))
)

microbenchmark::microbenchmark(
  ldmvnorm = ldmvnorm(t(obs), chol = C, logLik = FALSE),
  dmvnorm = apply(obs, 1, function(x) dmvnorm(x, sigma = S, log = TRUE)),
  times = 10
) |> plot()

ldmvnorm(t(obs) + m, mean = m, chol = C, logLik = TRUE) |> str()
sldmvnorm(t(obs) + m, mean = m, chol = C, logLik = TRUE) |> lapply(unclass) |> str()

###### CALCOLO DELLA VEROSIMIGLIANZA #######

J <- 5
set.seed(2121)
R <- sample(c(outer(c(-1, 0, 0, .35, 1), c(.25, .5, .75))), J^2, replace = TRUE) |> matrix(nrow = J)
R[upper.tri(R)] <- t(R)[upper.tri(R)]
diag(R) <- 1
Sigma <- diag(sqrt(1:J / 2)) %*% R %*% diag(sqrt(1:J / 2))
C <- t(chol(Sigma))

# conversione in ltMatrices

prm <- C[lower.tri(C, diag = TRUE)]
lt <- ltMatrices(matrix(prm, ncol = 1L),
                   diag = TRUE, ### has diagonal elements
                   byrow = FALSE) ### prm is column-major
all.equal(C, as.array(lt)[,,1], check.attributes = FALSE, check.class = FALSE)
all.equal(Sigma, as.array(Tcrossprod(lt))[,,1], check.attributes = FALSE, check.class = FALSE)

# Vogliamo stimare C in max verosimiglianza
# We generate some data from N (0 , Σ) by first sampling from Z ~ N (0J , I ) and then
# computing Y = CZ + μ ∼ NJ (μ, CC^⊤)
N <- 200
mn <- 1:J # medie
Z <- matrix(rnorm(N * J), nrow = J)
Y <- Mult(lt, Z) + mn

# stime di massima verosimiglianza analitiche
(muhat <- rowMeans(Y))
(Shat <- var(t(Y)) * (N - 1) / N)

# stime numeriche
Yc <- Y - rowMeans(Y) # centratura della Y
ll <- function(parm) { # -loglik
  Cc <- ltMatrices(parm, diag = TRUE, byrow = FALSE)
  -ldmvnorm(obs = Yc, chol = Cc)
}
sc <- function(parm) { # score function
  Cc <- ltMatrices(parm, diag = TRUE, byrow = FALSE)
  -rowSums(unclass(sldmvnorm(obs = Yc, chol = Cc)$chol))
}
# verifichiamo che llik e score funzionino
cML <- t(chol(Shat))[lower.tri(Shat, diag = TRUE)]
ll(cML)
start <- runif(length(cML))
all.equal(numDeriv::grad(ll, start), sc(start), check.attributes = FALSE, check.class = FALSE)

# ottimizzazione
# vincoli per i parametri di varianza (sulla diagonale)
llim <- rep(-Inf, J * (J + 1) / 2)
llim[which(rownames(unclass(lt)) %in% paste(1:J, 1:J, sep = "."))] <- 1e-4
op <- optim(start, fn = ll, gr = sc, method = "L-BFGS-B",
            lower = llim, control = list(trace = FALSE))
all.equal(op$par, cML) ## OK

# aggiungiamo anche la stima MV di mu

ll2 <- function(parm) { # -loglik
  muc <- parm[1:J]
  Cc <- ltMatrices(parm[-(1:J)], diag = TRUE, byrow = FALSE)
  -ldmvnorm(obs = Y, mean = muc, chol = Cc)
}
sc2 <- function(parm) { # score function
  muc <- parm[1:J]
  Cc <- ltMatrices(parm[-(1:J)], diag = TRUE, byrow = FALSE)
  -c(Mult(Crossprod(solve(Cc)), rowSums(Y - muc)), # score analitica per la media gaussiana V^(-1) sum(y_j - mu)
     rowSums(unclass(sldmvnorm(obs = Y, mean = muc, chol = Cc)$chol)))
}
cML2 <- c(unname(muhat), cML)
ll2(cML2)
sc2(cML2)
start2 <- runif(length(cML2))
all.equal(grad(ll2, start2), sc2(start2), check.attributes = FALSE, check.class = FALSE)
llim2 <- c(rep(-Inf, J), llim)
op2 <- optim(start2, fn = ll2, gr = sc2, method = "L-BFGS-B",
            lower = llim2, control = list(trace = FALSE, factr = 1))
all.equal(op2$par, cML2) ## OK
c(ll2(op2$par), ll2(cML2))

# riparametrizzazione con la varianza

ll3 <- function(parm) { # -loglik
  muc <- parm[1:J]
  Cc <- try(chol(ss <- syMatrices(parm[-(1:J)], diag = TRUE, byrow = FALSE)),
            silent = TRUE)
  if (inherits(Cc, "try-error"))
    return(1e10)
  -ldmvnorm(obs = Y, mean = muc, chol = Cc)
}
sc3 <- function(parm) { # score function
  muc <- parm[1:J]
  S <- as.array(syMatrices(parm[-(1:J)], diag = TRUE, byrow = FALSE))[, , 1]
  Scinv <- solve(S)
  # score analitica per la media gaussiana V^(-1) sum(y_j - mu)
  scmu <- Scinv %*% rowSums(Y - muc)
  # score analitica per la varianza gaussiana
  n <- ncol(Y)
  scQ <- vapply(seq_len(ncol(Y)),
                \(i) Scinv %*% tcrossprod(Y[, i] - muc) %*% Scinv,
                FUN.VALUE = array(0, dim(Scinv)))
  scQ <- apply(scQ, 1:2, sum)
  scoreS <- -n * Scinv + n/2 * diag(diag(Scinv)) + scQ - 0.5 * diag(diag(scQ))
  -c(scmu, scoreS[lower.tri(scoreS, diag = TRUE)])
}
# controllo anche la versione implementata da me nel pacchetto

sc32 <- function(parm) {
  muc <- parm[1:J]
  S <- as.array(syMatrices(parm[-(1:J)], diag = TRUE, byrow = FALSE))[, , 1]
  sc <- apply(Y, 2, .score1mat, mu = muc, Sigma = S)
  -rowSums(vapply(sc,
                  \(x) c(x[[1]], x[[2]][lower.tri(x[[2]], diag = TRUE)]),
                  numeric(length(parm))))
}
cML3 <- c(unname(muhat), Shat[lower.tri(Shat, diag = TRUE)])
ll3(cML3)
sc3(cML3)
all.equal(sc3(cML3), sc32(cML3))
start3 <- c(rep(0, J), diag(J)[lower.tri(diag(J), diag = TRUE)])
all.equal(numDeriv::grad(ll3, start3), sc3(start3), check.attributes = FALSE, check.class = FALSE)
llim3 <- c(rep(-Inf, J), llim)
op3 <- optim(start3, fn = ll3, gr = sc3, method = "L-BFGS-B",
             lower = llim3, control = list(trace = FALSE))
all.equal(op3$par, cML3) ## OK
c(ll3(op3$par), ll3(cML3)) |> as.character()

# trasformazione logaritmica dei parametri di varianza di una normale

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
cML4 <- cML3
cML4[J + whichlog] <- log(cML4[J + whichlog])
all.equal(ll4(cML4), ll3(cML3))
all.equal(sc4(cML4), sc3(cML3))
start4 <- start3
start4[J + whichlog] <- log(start4[J + whichlog])
all.equal(ll4(start4), ll3(start3))
all.equal(sc4(start4), sc3(start3))

(op4 <- nlminb(start4, ll4, sc4))
(op42 <- nlminb(start4, ll4))
 op4par <- op4$par
op4par[J + whichlog] <- exp(op4par[J + whichlog])
all.equal(op4par, op3$par)


# trasformazione logaritmica dei parametri di varianza, con trasformazione dei parametri

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
cML5 <- c(mu0 = cML3[1], b = cML3[2:J] / cML3[1],
          s = optim(rep(1, J),
                    \(x) sum((Shat -
                                diagoffdiag(1, .5, J) * tcrossprod(sqrt(exp(x))))^2),
                    control = list(maxit = 2e6))$par)
cMLbis <- c(cML3[1:J], (diagoffdiag(1, .5, J) * tcrossprod(sqrt(exp(cML5[-(1:J)]))))[lower.tri(diag(J), diag = TRUE)])
all.equal(ll5(cML5), ll3(cMLbis))
start5 <- rnorm(length(cML5))

(op5 <- nlminb(start5, ll5, sc5))
(op52 <- nlminb(start5, ll5, control = list(eval.max = 30000, iter.max = 30000)))
op4par <- op4$par
op4par[J + whichlog] <- exp(op4par[J + whichlog])
all.equal(op4par, op3$par)



microbenchmark::microbenchmark(
  grad = nlminb(start4, ll4, sc4),
  numderiv = nlminb(start4, ll4, \(x) numDeriv::grad(ll4, x)),
  default = nlminb(start4, ll4),
  times = 10
) |> plot()
