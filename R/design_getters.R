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
