.build.llik <- function(llik.fun, score.fun = NULL, use.data, stop.on.fail) {
  GETDATA <- list(
    y = substitute(crr.get.theta(DATA, raw = TRUE),
                   list(DATA = if (use.data) quote(data) else quote(object))),
    Gamma = substitute(crr.get.Gamma(DATA, raw = TRUE),
                       list(DATA = if (use.data) quote(data) else quote(object)))
  )
  TERMINATE_NA <-if (is.na(stop.on.fail)) quote({
    return(-Inf)
  }) else if (stop.on.fail) quote({
    stop("rilevati NA nella verosimiglianza")
  }) else quote({
    warning("Si sono prodotti NA nel calcolo della verosimiglianza; ",
            "si restituisce -Inf.")
    return(-Inf)
  })
  pieces <- list(
    `__GETPARS__` = quote({
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
      dun.pars <- lapply(dun, \(d) lapply(params, \(x) if (length(x) == 1) x else x[d]))
      mu <- .mapply(crr.mu.int, list(dun.pars, dbs), NULL)[dun.map]
      Sigma <- .mapply(crr.sigma.int, list(dun.pars, dbs), NULL)[dun.map]
    }),
    `__CHECK_AND_LOG__` = substitute({
      if (anyNA(ll)) TERMINATE_NA
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
    }, list(TERMINATE_NA = TERMINATE_NA))
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
  body(llik.fun) <- eval(substitute(substitute(EXPR, pieces),
                                    list(EXPR = body(llik.fun))))
  if (!is.null(score.fun)) {
    pieces$`__LLIKFN__` <- substitute(llik <- FUN, list(FUN = llik.fun))
    formals(score.fun) <- arglist
    body(score.fun) <- eval(substitute(substitute(EXPR, pieces),
                                       list(EXPR = body(score.fun))))
  }
  if (is.null(score.fun))
    return(llik.fun)
  structure(llik.fun, score = score.fun)
}

##'
##' A partire da un *design*, crea un oggetto funzione che può essere
##' passato ad `optim`.
##'
##' @rdname llik-from-design
##' @export
##' @title Funzione per la definizione della verosimiglianza per NCRR.
##'
##' @param object oggetto che definisce il design
##' @param transform applica trasformazioni ai parametri di
##'   varianza/correlazione
##' @param echo regola il livello delle stampe di debug
##' @param vcov.type struttura della matrice di varianza-covarianza
##' @param use.data se `TRUE`, la funzione risultante prenderà in input il dataset
##'   completo anzichè le singole componenti. Utile per l'utilizzo in combinazione
##'   con la libreria `boot` o `likelihoodAsy`.
##' @param stop.on.fail La funzione dovrebbe restituire un errore se una
##'   componente della verosimiglianza risulta NA? Se l'argomento è impostato a
##'   FALSE oppure NA, la funzione restituirà invece -Inf, nel primo caso con un
##'   avvertimento.
##' @return una funzione del vettore dei parametri (in tal senso `llik1` ne
##'   costituisce una versione semplificata), con la possibilità di fissare gli
##'   stessi (opzione `fixed = list(...)`)
##' @author Marco Bressan
get.llik.from.design <- function(object, transform = TRUE, echo = 0,
                                 vcov.type = attr(object, "vcov.type"),
                                 use.data = FALSE, stop.on.fail = FALSE) {
  np <- length(tt <- unique(do.call(c, object$design))) - 1
  if (is.null(vcov.type))
    vcov.type <- "normal"
  vcov.type <- match.vcov.type(vcov.type)
  fixed.default <- match.vcov.fixed(vcov.type, TRUE, np)
  # quantità precalcolate
  n <- length(object$design)
  dun <- unique(object$design)
  dun.map <- match(object$design, dun)
  dbs <- vapply(dun, \(d) 0 %in% d, logical(1))
  par.pos <- crr.par.idx(np, fixed = names(fixed.default),
                         parlen = attr(fixed.default, "parlen"))
  par.trans <- grep("sigma|rho", names(par.pos))
  heps <- .Machine$double.eps^(1/3)
  # costruzione likelihood
  .build.llik(
    llik.fun = function(params, y = crr.get.theta(object, raw = TRUE),
                        Gamma = crr.get.Gamma(object, raw = TRUE),
                        fixed = NULL) {
      pold <- params
      `__GETPARS__`
      ll <- mapply(\(t, m, Si, Gi) {
        return(mvtnorm::dmvnorm(t, m, Si + Gi, log = TRUE))
        # tengo la parte seguente solo per informazione
        ## chl <- try(chol(S <- Si + Gi), silent = echo <= 3)

        ## if (inherits(chl, "try-error")) {
        ##   #browser()
        ##   chl <- S * NaN
        ##   #chl <- as.matrix(as(Matrix::Cholesky(S),"dtrMatrix"))
        ## }
        ## if (!exists(".__likdbg", globalenv()))
        ##   assign(".__likdbg", list(list(S, chl)), globalenv())
        ## else
        ##   assign(".__likdbg",
        ##          append(get(".__likdbg", globalenv()), list(list(S, chl))), globalenv())
        ## chl <- mvtnorm::ltMatrices(chl[which(upper.tri(chl, diag = TRUE))], diag = TRUE)
        ## ldn <- mvtnorm::ldmvnorm(t, mean = m, chol = chl)
        ## mvn <- mvtnorm::dmvnorm(t, m, S, log = TRUE)
        ## mvn
      }, y, mu, Sigma, Gamma)
      `__CHECK_AND_LOG__`
      return(ll)
    },
    score.fun <- function(params, y = crr.get.theta(object, raw = TRUE),
                          Gamma = crr.get.Gamma(object, raw = TRUE),
                          fixed = NULL) {
      gr <- params * 0
      pold <- params
      #`__LLIKFN__` # llik <- function(...){...}
      ll <- matrix(NA, 2, n)
      for (i in seq_along(gr)) {
        for (l in c(-1, 1)) {
          params <- pold
          params[i] <- params[i] + l * heps
          `__GETPARS__`
          ll[as.integer(l / 2 + 1.5), ] <- mapply(\(t, m, Si, Gi) {
            return(mvtnorm::dmvnorm(t, m, Si + Gi, log = TRUE))
          }, y, mu, Sigma, Gamma)
        }
        gr[i] <- sum(diff(ll) / 2 / heps)
      }
      gr
    },
    use.data = as.logical(use.data),
    stop.on.fail = as.logical(stop.on.fail)
  )
}

##' @rdname llik-from-design
##' @details
##' `get.llik.from.design2` è una specializzazione per design contenenti solo studi
##' con due trattamenti.
##'
get.llik.from.design2 <- function(object, transform = TRUE, echo = 0,
                                 vcov.type = attr(object, "vcov.type"),
                                 use.data = FALSE, stop.on.fail = TRUE) {
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
        object = mapply(\(Si, Gi) {
          chl <- try(chol(S <- Si + Gi), silent = echo <= 3)
          if (inherits(chl, "try-error"))
            chl <- as.matrix(as(Matrix::Cholesky(S),"dtrMatrix"))
          chl[which(upper.tri(Si, diag = TRUE))]
        } , Sigma, Gamma),
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
    use.data = as.logical(use.data),
    stop.on.fail = as.logical(stop.on.fail)
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

#' @export
vcov.ncrr.design <- function(object, x0, llik.fn = get.llik.from.design(object),
                             score.fn = attr(llik.fn, "score"), sandwich = FALSE,
                             ...,
                             score = NULL, J = NULL) {
  n <- length(object$design)
  J <- J %||% -optimHess(x0, llik.fn, gr = score.fn, ...)
  invJ <- solve(J)
  if (isFALSE(sandwich))
    return(invJ)
  score <- if (is.null(score.fn)) pracma::grad(llik.fn, x0) else score.fn(x0)
  I <- tcrossprod(score)
  invJ %*% I %*% invJ
}
