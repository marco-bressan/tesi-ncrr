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

    ll <- mapply(\(t, m, s, g) {
      mvtnorm::dmvnorm(t, mean = m, sigma = s + g, log = TRUE)
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
