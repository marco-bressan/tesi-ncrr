### five observations
N <- 5000L
### dimension
J <- 4


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
all.equal(grad(ll, start), sc(start), check.attributes = FALSE, check.class = FALSE)

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
  Cc <- try(chol(ss <- syMatrices(parm[-(1:J)], diag = TRUE, byrow = FALSE)))
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
all.equal(grad(ll3, start3), sc3(start3), check.attributes = FALSE, check.class = FALSE)
llim3 <- c(rep(-Inf, J), llim)
op3 <- optim(start3, fn = ll3, gr = sc3, method = "L-BFGS-B",
             lower = llim3, control = list(trace = FALSE))
all.equal(op3$par, cML3) ## OK
c(ll3(op3$par), ll3(cML3)) |> as.character()

