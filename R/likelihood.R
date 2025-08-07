##' Funzione per l'ottimizzazione centralizzata
##'
##' @title Ottimizzazione di una network crr.
##' @param design Disegno dello studio (creato con l'apposita funzione)
##' @param init Valori iniziali in formato lista. I parametri esplicitamente
##'   settati a `NULL` sono assunti fissati; quelli mancanti vengono inferiti
##'   automaticamente.
##' @param rs Numero di random-start. Se impostato, sovrascrive i valori di
##'   partenza specificati da init
##' @return Risultato dell'ottimizzazione
##' @author Marco Bressan
crr.optim <- function(design, init, rs = 1, transform = TRUE,
                      optim.method = "Nelder-Mead",
                      optim.control = list(), ...) {
  .NotYetImplemented()
  fn <- get.llik.from.design(design, ...)
  stopifnot("I parametri `init` non sono nominati!" = !is.null(names(init)))
  #match(names(init),PNAMES) #???
  fixpars <- names(Filter(is.null, init))
  explicpars <- Filter(is.numeric, init)
  par.init <- getInitial(design, rep = rs, fixed = fixpars)
  opt.fn <- get.llik.from.design(design, transform = transform, )
  if (length(explicpars) > 0) {
    if (rs == 1) {
      par.init[names(explicpars)] <- explicpars
      oo <- optim(par.init, \(x) -opt.fn(x), control = optim.control)
    } else {
      mm <- outer(rep(1, rs), unlist(explicpars))
      ide <- grep(colnames(par.init), paste(names(explicpars), collapse = "|"))
      stopifnot(ncol(mm) == length(ide))
      par.init[, ide] <- mm
      #...
    }
  }
  #TODO: DA FINIRE
  #....
}

crr.fix.pars <- function(..., parname, value, np) {
  .NotYetImplemented()
  # da rifinire...
  parnames <- c("alpha", "beta", "mu0", "sigma20", "rho", "sigma2")
  parregex <- paste(parnames, collapse = "|")
  nm <- stringr::str_extract_all(sprintf("(%s)([0-9]+)", parregex), ...names())
  browser()
  params
}

##' @rdname llik-from-design
##' A partire da un *design*, crea un oggetto funzione che può essere
##' passato ad `optim`.
##'
##' @export
##' @title Funzione per la definizione della verosimiglianza per NCRR.
##' @param object oggetto che definisce il design
##' @param transform applica trasformazioni ai parametri di
##'   varianza/correlazione
##' @param echo regola il livello delle stampe di debug
##' @param vcov.type struttura della matrice di varianza-covarianza.
##' @return una funzione del vettore dei parametri (in tal senso `llik1` ne
##'   costituisce una versione semplificata), con la possibilità di fissare gli
##'   stessi (opzione `fixed = list(...)`)
##' @author Marco Bressan
get.llik.from.design <- function(object, transform = TRUE, echo = 0,
                                 vcov.type = attr(object, "vcov.type")) {
  np <- length(tt <- unique(do.call(c, object$design)))
  fixed.default <- NULL
  if (is.null(vcov.type))
    vcov.type <- "normal"
  vcov.type <- match.vcov.type(vcov.type)
  fixed.default <- match.vcov.fixed(vcov.type, TRUE, np - 1)
  GETPARS <- function(params, fixed) {
    if (length(names(fixed.default)) > 0) {
      fixed[names(fixed.default)] <- fixed.default
      attributes(fixed) <- attributes(fixed.default)
    }
    params <- crr.split.par(params, np - 1, transform = transform,
                            fixed = names(fixed),
                            parlen = attr(fixed, "parlen"))
    #params <- c(params, fixed)
    #params <- params[match(PNAMES, names(params))]
    params
  }

  # questo è il return,  la funzione obiettivo
  function(params, y = crr.get.theta(object, raw = TRUE),
           Gamma = crr.get.Gamma(object, raw = TRUE),
           fixed = NULL) {
    params <- GETPARS(params, fixed)
    if (echo > 1) {
      mapply( \(x, nm) paste(nm, "=", deparse1(round(x, 6))),
             params, names(params)) |>
        paste(collapse = ", ") |>
        cat("\n")
    }
    mu <- do.call(crr.get.mu,
                  append(list(object = object, raw = TRUE),
                         params[c("alpha", "beta", "mu0")]))
    Sigma <- do.call(crr.get.sigma,
                     append(list(object = object, raw = TRUE, type = vcov.type),
                            params[c("beta", "sigma20", "sigma2", "rho")]))
    ll <- mapply(\(t, m, Si, Gi) {
      Ci <- chol(Si + Gi)
      Ci <- mvtnorm::ltMatrices(Ci[which(upper.tri(Ci, diag = TRUE))], diag = TRUE)
      mvtnorm::ldmvnorm(t, mean = m, chol = Ci)
    }, y, mu, Sigma, Gamma)
    if (anyNA(ll)) {
      warning("Si sono prodotti NA nel calcolo della verosimiglianza,",
              " che sono stati scartati.")
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
    #browser()
    ll <- sum(ll, na.rm = TRUE)
    if (echo > 0)
      cat("CURRENT VALUE: ", ll, "\n")
    if (echo > 1)
      cat("=====================================\n")
    return(ll)
  }
}

##' @rdname llik-from-design
##' @details
##' `get.llik.from.design2` è una specializzazione per design contenenti solo studi
##' con due trattamenti.
##'
get.llik.from.design2 <- function(object, transform = TRUE, echo = 0,
                                 vcov.type = attr(object, "vcov.type")) {
  stopifnot("Tutti gli studi devono confrontare esattamente due trattamenti!" =
              lengths(object$design) == 2)
  np <- length(tt <- unique(do.call(c, object$design)))
  fixed.default <- NULL
  if (is.null(vcov.type))
    vcov.type <- "normal"
  vcov.type <- match.vcov.type(vcov.type)
  fixed.default <- match.vcov.fixed(vcov.type, TRUE, np - 1)
  GETPARS <- function(params, fixed) {
    if (length(names(fixed.default)) > 0) {
      fixed[names(fixed.default)] <- fixed.default
      attributes(fixed) <- attributes(fixed.default)
    }
    params <- crr.split.par(params, np - 1, transform = transform,
                            fixed = names(fixed),
                            parlen = attr(fixed, "parlen"))
    #params <- c(params, fixed)
    #params <- params[match(PNAMES, names(params))]
    params
  }

  # questo è il return, la funzione obiettivo e la sua score come attributo
  structure(
    function(params, y = crr.get.theta(object, raw = TRUE),
             Gamma = crr.get.Gamma(object, raw = TRUE),
             fixed = NULL) {
      params <- GETPARS(params, fixed)
      if (echo > 1) {
        mapply( \(x, nm) paste(nm, "=", deparse1(round(x, 6))),
                params, names(params)) |>
          paste(collapse = ", ") |>
          cat("\n")
      }
      mu <- do.call(crr.get.mu,
                    append(list(object = object, raw = TRUE),
                           params[c("alpha", "beta", "mu0")]))
      Sigma <- do.call(crr.get.sigma,
                       append(list(object = object, raw = TRUE, type = vcov.type),
                              params[c("beta", "sigma20", "sigma2", "rho")]))
      cholSigt <- mvtnorm::ltMatrices(
        object = mapply(\(Si, Gi) chol(Si + Gi)[which(upper.tri(Si, diag = TRUE))], Sigma, Gamma),
        diag = TRUE)
      tt <- do.call(cbind, y)
      mu <- do.call(cbind, mu)
      ll <- mvtnorm::ldmvnorm(tt, mean = mu, chol = cholSigt, logLik = FALSE)
      if (anyNA(ll)) {
        warning("Si sono prodotti NA nel calcolo della verosimiglianza,",
                " che sono stati scartati.")
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
      #browser()
      ll <- sum(ll, na.rm = TRUE)
      if (echo > 0)
        cat("CURRENT VALUE: ", ll, "\n")
      if (echo > 1)
        cat("=====================================\n")
      return(ll)
    },

    score = function(params, y = crr.get.theta(object, raw = TRUE),
                     Gamma = crr.get.Gamma(object, raw = TRUE),
                     fixed = NULL) {
      params <- GETPARS(params, fixed)
      if (echo > 1) {
        mapply( \(x, nm) paste(nm, "=", deparse1(round(x, 6))),
                params, names(params)) |>
          paste(collapse = ", ") |>
          cat("\n")
      }
      mu <- do.call(crr.get.mu,
                    append(list(object = object, raw = TRUE),
                           params[c("alpha", "beta", "mu0")]))
      Sigma <- do.call(crr.get.sigma,
                       append(list(object = object, raw = TRUE, type = vcov.type),
                              params[c("beta", "sigma20", "sigma2", "rho")]))
      # score rispetto ai paramentri normali (mu, Sigmatilde)
      scL <- mapply(.score1mat, y, mu, Sigma, SIMPLIFY = FALSE)
      # derivata di mu in funzione di alpha, beta, mu0
      scmu <- do.call(crr.get.scmu,
                      append(list(object = object),
                             params[c("beta", "mu0")]))
      # derivata di Sigmatilde in funzione di beta, sigma20, rho, sigma2
      scSigt <- do.call(crr.get.scSigma,
                        append(list(object = object),
                               params[c("beta", "sigma20", "rho", "sigma2")]))
      # regola del concatenamento (funziona anche per le matrici?!)
      sc <- mapply(\(scLi, scmui, scSigti) {
        c(alpha = crossprod(scLi$mu, scmui$alpha),
          beta = crossprod(scLi$mu, scmui$beta),
          mu0 = crossprod(scLi$mu, scmui$mu0),
          sigma20 = crossprod(c(scLi$Sigma), c(scSigti$sigma20)),
          rho = crossprod(c(scLi$Sigma), c(scSigti$rho)),
          sigma2 = crossprod(c(scLi$Sigma), apply(scSigti$sigma2, 3, c)))
          # dimensioni:    n_i x (k_i)^2  ,   (k_i)^2 x dim(sigma2)
      }, scL, scmu, scSigt)
      return(rowSums(sc))
    }
  )
}

.score1mat <- function(tt, mu, Sigma) {
  Scinv <- solve(Sigma)
  # score analitica per la media gaussiana V^(-1) sum(y_j - mu)
  scmu <- Scinv %*% (tt - mu)
  # score analitica per la varianza gaussiana
  scQ <- Scinv %*% tcrossprod(tt - mu) %*% Scinv
  scoreS <- -Scinv + 0.5 * diag(diag(Scinv)) + scQ - 0.5 * diag(diag(scQ))
  list(mu = scmu, Sigma = scoreS)
}


##' Verosimiglianza per un singolo studio per NCRR __con baseline__.
##' --VERSIONE PRELIMINARE--
##'
##' @title Verosimiglianza NCRR
##' @param params vettore dei parametri contenuti in un unico vettore; in
##'   ordine: alfa_0j, beta_0j, mu0, sigma^2_0, rho, sigma^2_{0j}.
##' @param y vettore dei vartheta
##' @param Gamma matrice Gamma
##' @param design vettore che specifica il design (attualmente utile solo la
##'   sua dimensione)
##' @return valore della log-verosimiglianza
##' @author Marco Bressan
##' @export
llik1 <- function(params, y, Gamma, design = c(0, 1)) {
  stop("Non usare questa funzione!")
  params <- crr.split.par(params, np <- length(design) - 1)
  mu <- rep(NA, np + 1)
  Sigma <- matrix(0, np + 1, np + 1)
  ll <- with(params, {
    mvtnorm::dmvnorm(
      y, mean = mu <<- crr.mean.baseline0(alpha, beta, mu0, design),
      sigma = Sigma <<- crr.vcov(beta, sigma20, sigma2, rho, design) + Gamma,
      log = TRUE)
  })
  #browser()
  cat("PARAMS: ")
  with(params, {
    lapply(list(alpha, beta, mu0, sigma20, rho, sigma2), \(x) deparse1(round(x, 6))) |>
      append(x = list(fmt = "alpha = %s, beta = %s, \n\tmu0 = %s, sigma20 = %s, rho = %s, \n\tsigma2 = %s"), values = _) |>
      do.call(sprintf, args = _) |>
      cat("\nVALUE: ", ll, "\n\n")
    print(list(marginal_mu = mu, marginal_tildeSigma = Sigma))
    cat("\n=====================================\n")
  })
  return(ll)
}
