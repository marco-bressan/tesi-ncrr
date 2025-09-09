.build.llik <- function(llik.fun, score.fun = NULL, use.data) {
  GETDATA <- list(
    y = substitute(crr.get.theta(DATA, raw = TRUE),
                   list(DATA = if (use.data) quote(data) else quote(object))),
    Gamma = substitute(crr.get.Gamma(DATA, raw = TRUE),
                       list(DATA = if (use.data) quote(data) else quote(object)))
  )
  pieces <- alist(
    `__GETPARS__` = {
      changed <- !is.null(fixed)
      if (length(names(fixed.default)) > 0) {
        fixed[names(fixed.default)] <- fixed.default
        attributes(fixed) <- attributes(fixed.default)
      }
      if (changed) {
        par.pos <- crr.par.idx(np, fixed = names(fixed), parlen = attr(fixed, "parlen"))
        par.trans <- grep("sigma|rho", names(par.pos))
      }
      params <- .parsplit1(params, par.pos, par.trans)
      if (length(fixed) > 0)
        params <- append(fixed, params)
      if (vcov.type != "normal")
        params <- set.vcov.params(params, np, vcov.type)
      if (echo > 1) {
        mapply( \(x, nm) paste(nm, "=", deparse1(round(x, 6))),
               params, names(params)) |>
          paste(collapse = ", ") |>
          cat("\n")
      }
      mu <- crr.get.mu(object, params, raw = TRUE)
      Sigma <- crr.get.sigma(object, params, raw = TRUE)
    },
    `__CHECK_AND_LOG__` = {
      if (anyNA(ll)) {
        warning("Si sono prodotti NA nel calcolo della verosimiglianza, ",
                "che sono stati scartati.")
        stop("rilevati NA nella verosimiglianza")
      }
      if (echo > 2) {
        cat("CURRENT PIECEWISE LLIK:\n")
        mapply(\(m, s, l) {
          colnames(s) <- c("SIGMA", rep("", ncol(s) - 1))
          print(cbind("MU" = m, s))
          print(c("LLIK" = l))
        }, mu, Sigma, ll)
      }
      ll <- sum(ll, na.rm = TRUE)
      if (echo > 0)
        cat("CURRENT VALUE: ", ll, "\n")
      if (echo > 1)
        cat("=====================================\n")
    }
  )
  arglist <- append(alist(params =, fixed = NULL),
                   if (use.data) alist(data = object) else GETDATA,
                   after = 1)
  if (use.data) {
    GETDATAexpr <- append(as.symbol("{"),
                          .mapply(\(expr, vname) call("<-", as.name(vname), expr),
                                  list(GETDATA, names(GETDATA)), NULL),
                          after = 1)
    pieces$`__GETPARS__` <- .append.expr(pieces$`__GETPARS__`, as.call(GETDATAexpr))
  }
  formals(llik.fun) <- arglist
  body(llik.fun) <- .subst.keys(body(llik.fun), pieces)
  if (!is.null(score.fun)) {
    formals(score.fun) <- arglist
    body(score.fun) <- .subst.keys(body(score.fun), pieces)
  }
  if (is.null(score.fun))
    return(llik.fun)
  structure(llik.fun, score = score.fun)
}

##' @rdname llik-from-design
##' A partire da un *design*, crea un oggetto funzione che può essere
##' passato ad `optim`.
##'
##' @export
##' @title Funzione per la definizione della verosimiglianza per NCRR.
##'
##' @param object oggetto che definisce il design
##' @param transform applica trasformazioni ai parametri di
##'   varianza/correlazione
##' @param echo regola il livello delle stampe di debug
##' @param use.data se `TRUE`, la funzione risultante prenderà in input il dataset
##' completo anzichè le singole componenti. Utile per l'utilizzo in combinazione
##' con la libreria `boot` o `likelihoodAsy`.
##' @param vcov.type struttura della matrice di varianza-covarianza.
##'
##' @return una funzione del vettore dei parametri (in tal senso `llik1` ne
##'   costituisce una versione semplificata), con la possibilità di fissare gli
##'   stessi (opzione `fixed = list(...)`)
##' @author Marco Bressan
get.llik.from.design <- function(object, transform = TRUE, echo = 0,
                                 vcov.type = attr(object, "vcov.type"),
                                 use.data = FALSE) {
  np <- length(tt <- unique(do.call(c, object$design))) - 1
  fixed.default <- NULL
  if (is.null(vcov.type))
    vcov.type <- "normal"
  vcov.type <- match.vcov.type(vcov.type)
  fixed.default <- match.vcov.fixed(vcov.type, TRUE, np)
  par.pos <- crr.par.idx(np, fixed = names(fixed.default), parlen = attr(fixed.default, "parlen"))
  par.trans <- grep("sigma|rho", names(par.pos))
  .build.llik(
    function(params, y = crr.get.theta(object, raw = TRUE),
             Gamma = crr.get.Gamma(object, raw = TRUE),
             fixed = NULL) {
      `__GETPARS__`
      ll <- mapply(\(t, m, Si, Gi) {
        Ci <- chol(Si + Gi)
        Ci <- mvtnorm::ltMatrices(Ci[which(upper.tri(Ci, diag = TRUE))], diag = TRUE)
        mvtnorm::ldmvnorm(t, mean = m, chol = Ci)
      }, y, mu, Sigma, Gamma)
      `__CHECK_AND_LOG__`
      return(ll)
    },
    use.data = as.logical(use.data)
  )
}

##' @rdname llik-from-design
##' @details
##' `get.llik.from.design2` è una specializzazione per design contenenti solo studi
##' con due trattamenti.
##'
get.llik.from.design2 <- function(object, transform = TRUE, echo = 0,
                                 vcov.type = attr(object, "vcov.type"),
                                 use.data = FALSE) {
  stopifnot("Tutti gli studi devono confrontare esattamente due trattamenti!" =
              lengths(object$design) == 2)
  np <- length(tt <- unique(do.call(c, object$design))) - 1
  fixed.default <- NULL
  if (is.null(vcov.type))
    vcov.type <- "normal"
  vcov.type <- match.vcov.type(vcov.type)
  fixed.default <- match.vcov.fixed(vcov.type, TRUE, np)
  sctrans <- list(sigma = \(x) exp(x),
                  rho = \(x) 4 * exp(2*x) / (exp(2*x) + 1)^2)
  par.pos <- crr.par.idx(np, fixed = names(fixed.default), parlen = attr(fixed.default, "parlen"))
  par.trans <- grep("sigma|rho", names(par.pos))
  # questo è il return, la funzione obiettivo e la sua score come attributo
  .build.llik(
    llik.fun = function(params, y = crr.get.theta(object, raw = TRUE),
             Gamma = crr.get.Gamma(object, raw = TRUE),
             fixed = NULL) {
      `__GETPARS__`
      cholSigt <- mvtnorm::ltMatrices(
        object = mapply(\(Si, Gi) chol(Si + Gi)[which(upper.tri(Si, diag = TRUE))], Sigma, Gamma),
        diag = TRUE)
      tt <- do.call(cbind, y)
      mu <- do.call(cbind, mu)
      ll <- mvtnorm::ldmvnorm(tt, mean = mu, chol = cholSigt, logLik = FALSE)
      `__CHECK_AND_LOG__`
      return(ll)
    },
    score.fun =  function(params, y = crr.get.theta(object, raw = TRUE),
                     Gamma = crr.get.Gamma(object, raw = TRUE),
                     fixed = NULL) {
      onames <- names(params)
      opars <- params
      `__GETPARS__`
      # score rispetto ai paramentri normali (mu, Sigmatilde)
      scL <- mapply(.score1mat, y, mu, Sigma, SIMPLIFY = FALSE)
      # derivata di mu in funzione di alpha, beta, mu0
      scmu <- crr.get.scmu(object, params)
      # derivata di Sigmatilde in funzione di beta, sigma20, rho, sigma2
      scSigt <- crr.get.scSigma(object, params)
      # regola del concatenamento (funziona anche per le matrici?!)
      sc <- mapply(\(scLi, scmui, scSigti) {
        ret <- c(alpha = crossprod(scLi$mu, scmui$alpha),
          beta = crossprod(scLi$mu, scmui$beta),
          mu0 = crossprod(scLi$mu, scmui$mu0),
          sigma20 = crossprod(c(scLi$Sigma), c(scSigti$sigma20)),
          rho = crossprod(c(scLi$Sigma), c(scSigti$rho)),
          sigma2 = crossprod(c(scLi$Sigma), apply(scSigti$sigma2, 3, c)))
          # dimensioni:    n_i x (k_i)^2  ,   (k_i)^2 x dim(sigma2)
        if (isTRUE(transform)) { # derivata della trasformazione
          spars.ind <- grep("sigma", names(ret), value = TRUE)
          ret[spars.ind] <- ret[spars.ind] * sctrans$sigma(opars[spars.ind])
          ret["rho"] <- sctrans$rho(ret["rho"])
        }
        ret
      }, scL, scmu, scSigt)
      if (!is.null(onames))
        sc <- sc[onames, ] # toglie parametri inutili tipo rho per "achana"
      return(rowSums(sc))
    },
    use.data = as.logical(use.data)
  )
}

.score1mat <- function(tt, mu, Sigma) {
  Scinv <- solve(Sigma)
  scmu <- Scinv %*% (tt - mu)
  # score analitica per la varianza gaussiana
  scQ <- Scinv %*% tcrossprod(tt - mu) %*% Scinv
  scoreS <- -Scinv + 0.5 * diag(diag(Scinv)) + scQ - 0.5 * diag(diag(scQ))
  list(mu = scmu, Sigma = scoreS)
}
