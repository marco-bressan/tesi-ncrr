match.vcov.type <- function(type = c("normal", "achana", "equivar", "simple")) {
  match.arg(type, several.ok = FALSE)
}

set.vcov.params <- function(params, np = length(params$beta), vcov.type) {
  if (vcov.type == "achana") {
    params$sigma2 <- rep(params$sigma2[1], np)
    params$rho <- 0.5
  } else if (vcov.type == "equivar") {
    params$sigma2 <- rep(params$sigma20, np)
  } else if (vcov.type == "simple") {
    params$sigma2 <- rep(params$sigma20, np)
    params$rho <- 0.5
  }
  params
}

match.vcov.fixed <- function(type, value = FALSE, np) {
  type <- match.vcov.type(type)
  ff <- switch(type, achana = "rho", equivar = "sigma2", simple = c("rho", "sigma2"))
  if (isTRUE(value)) {
    ff <- as.list(setNames(nm = ff))
    if (type == "simple") {
      ff[["sigma2"]] <- rep(NA, np)
      ff[["rho"]] <- NA
    } else if (type == "achana") {
      ff[["rho"]] <- .5
    } else if (type == "equivar") {
      ff[["sigma2"]] <- rep(NA, np)
    }
  }
  if (type == "achana") {
    attr(ff, "parlen") <- c("sigma2" = 1)
  }
  ff
}

crr.get.sigma <- function(object, params, raw = FALSE) {
  dd <- object$design
  dunique <- unique(dd)
  #browser()
  Sigmal <- lapply(dunique, \(d) {
    psel <- par.select.multi(d, params)
    V <- crr.vcov(psel, d)
    if (anyNA(V)) stop("Na rilevati nel calcolo di sigma!") # togliere per efficientamento
    V
  })
  if (raw)
    return(Sigmal[match(dd, dunique)])
  blockdiag(Sigmal[match(dd, dunique)])
}

crr.get.mu <- function(object, params, raw = FALSE) {
  dd <- object$design
  #stopifnot("Baseline != 0 ancora da implementare!" = sapply(dd, \(d) d[1] == 0))
  dunique <- unique(dd)
  #browser()
  mul <- lapply(dunique, \(d) {
    #if (!0 %in% d) browser()
    psel <- par.select.multi(d, params)
    crr.mean(psel, d)
  })
  if (raw)
    return(mul[match(dd, dunique)])
  do.call(c, mul[match(dd, dunique)])
}

par.select.multi <- function(object, params) {
  lapply(params, \(x) if (length(x) == 1) x else x[object])
}

crr.get.Gamma <- function(object, raw = FALSE) {
  if (!raw)
    return(diag(object$gamma))
  lc <- cumsum(ll <- lengths(object$design))
  lc <- c(0, lc[-length(lc)])
  mapply(\(pos, len) diag(object$gamma[(pos + 1):(pos + len)]),
         lc, ll, SIMPLIFY = FALSE)
}

crr.get.theta <- function(object, raw = FALSE) {
  if (!raw)
    return(object$theta)
  lc <- cumsum(ll <- lengths(object$design))
  lc <- c(0, lc[-length(lc)])
  mapply(\(pos, len) object$theta[(pos + 1):(pos + len)],
         lc, ll, SIMPLIFY = FALSE)
}

##' L'uso della presente è confinato alle funzioni che calcolano i
##' parametri esplicitamente, contenute nella dir `data-raw/`
##'
##' @title Ottieni matrici dal design
##' @param object design di uno studio ncrr
##' @param what al momento solo "theta" o "gamma"
##' @return matrice di interesse (quella dei gamma o dei theta) in un
##'   formato che sia compatibile con la notazione usata in maxima nel calcolo
##'   delle derivate esplicite.
##' @author Marco Bressan
get.matrix.from.design <- function(object, what = c("theta", "gamma")) {
  what <- match.arg(what, several.ok = FALSE)
  stopifnot("passato oggetto non valido" = length(ll <- lengths(object$design)) > 0,
            "alcuni design sono nulli" = ll != 0)
  maxd <- max(unlist(object$design))
  mm <- matrix(NA, maxd + 1, length(object$design))
  for (i in seq_along(ll)) {
    mm[object$design[[i]] + 1, i] <- 1
  }
  stopifnot("lunghezza di `what` non compatibile" =
              sum(mm, na.rm = TRUE) == length(object[[what]]))
  mm[!is.na(mm)] <- object[[what]]
  t(mm)
}

crr.get.scmu <- function(object, params) {
  stopifnot(!USA_MIA_MODELLAZIONE)
  dd <- object$design
  #dd <- list(c(0, 1), c(0, 2), c(1, 3), c(0, 1, 2), c(2, 3))
  pos <- cumsum(lengths(dd))
  ll <- pos[length(pos)]
  faclev <- sort(unique(do.call(c, dd)))
  M0 <- numeric(ll)
  np <- length(params$beta)
  M1 <- matrix(0, ll, np)
  pos <- c(0, pos)
  lapply((dd), \(jj) {
    #jj <- dd[[i]]
    M <- vapply(jj, \(j) {
      row <- numeric(np)
      if (j > 0) row[j] <- 1
      row
    }, FUN.VALUE = numeric(np))
    if (all(jj > 0))
      M <- M[, c(2, 1)] - M # solo per due studi senza baseline
    M0 = as.numeric(jj == 0)
    M1 = t(M)
    list(alpha = M1, beta = params$mu0 * M1, mu0 = M0 + M1 %*% params$beta)
    #M0[(pos[i] + 1):pos[i + 1]] <- as.numeric(jj == 0)
    #M1[(pos[i] + 1):pos[i + 1], ] <- t(M)
  })
  #cbind(alpha = M1, beta = params$mu0 * M1, mu0 = M0 + M1 %*% params$beta)
}

crr.get.scSigma <- function(object, params) {
  dd <- object$design
  np <- length(params$beta)
  scs0 <- tcrossprod(c(1, params$beta))
  Dst <- diag(c(0, rep_len(sqrt(params$sigma2), np)))
  if (is.null(params$rho))  {stop("punto irraggiungibile!"); params$rho <- 0.5} #TODO: hardcoded
  Rst <- cbind(0, rbind(0, matrix(params$rho, np, np)))
  diag(Rst) <- c(0, rep_len(1, np))
  dRst <- Rst * 0
  dRst[-1, -1] <- 1 - diag(np)
  scs <- vapply(seq_along(params$sigma2), \(j) {
    dDstj <- Dst * 0
    dDstj[j + 1, j + 1] <- params$sigma2[j]
    (kronecker(Dst, dDstj) + kronecker(dDstj, Dst)) %*% c(Rst)
  }, FUN.VALUE = array(NA_real_, dim(Dst)))
  scrho <- Dst %*% dRst %*% Dst
  lapply (dd, \(d) {
    if (0 %in% d) {
      s0j <- scs0[d + 1, d + 1]
      rhoi <- scrho[d + 1, d + 1]
    } else {
      betai <- c(params$beta[rev(d)] - params$beta[d])
      s0j <- tcrossprod(betai)
      rhoi <- matrix(0, 2, 2) # con solo due studi sempre = 0 (teoricamente scrho[c(1, d[2]), c(1, d[2])])
    }
    list(sigma20 = s0j, rho = rhoi,
         sigma2 = scs[d + 1, d + 1 , if(length(params$sigma2 == 1))  1 else d,
                      drop = FALSE]) #TODO: ottimizzare; essenziale mantenere 3^ dim
  })
}
