blockdiag <- function(mats, fill = 0) {
  if (is.matrix(mats)) return(mats)
  stopifnot("only square matrices are supported!" = sapply(mats, nrow) == (nc <- sapply(mats, ncol)))
  out <- matrix(fill, sum(nc), sum(nc))
  ncc <- c(0, cumsum(nc))
  for (i in seq_along(mats)) {
    out[(ncc[i]+1):ncc[i + 1], (ncc[i]+1):ncc[i + 1]] <- mats[[i]]
  }
  out
}

diagoffdiag <- function(x = 1, offx = 0, n = length(x)) {
  offx * (1 - diag(n)) + x * diag(n)
}

##' Effettua una sostituzione nella lista dei parametri di una NCRR
##'
##'
##' @title Sostituisci parametri
##' @param params Lista o vettore dei parametri
##' @param subst Lista o vettore delle sostituzioni.
##' La sostituzione avviene in base al nome, che dunque è obbligatorio per
##' ogni parametro.
##' @return `params` opportunamente modificato
##' @author Marco Bressan
subst.params <- function(params, subst) {
  for (p in names(subst)) {
    if (length(params[[p]]) != length(subst[[p]])) {
      "Parametro %s fissato: sono stati forniti %i valori, ma ne sono richiesti %i" |>
        sprintf(p, length(subst[[p]]), length(params[[p]])) |>
        stop()
    }
    params[[p]] <- subst[[p]]
  }
  return(params)
}

.append.expr <- function(expr1, expr2, after = 1) {
  symb <- as.symbol("{")
  if (!identical(expr1[[1]], symb))
    expr1 <- substitute({EXPR}, list(EXPR = expr1))
  if (!identical(expr2[[1]], symb))
    expr2 <- substitute({EXPR}, list(EXPR = expr2))
  after <- min(max(after, 1), length(expr1))
  as.call(append(as.list(expr1), as.list(expr2)[-1], after = after))
}

.subst.keys <- function(expr, keys) {
  call.names <- vapply(expr, \(x) as.character(if (is.symbol(x)) x else x[[1]]),
                       character(1))
  subst.idx <- grep("^__.*__$", call.names)
  matched <- match(call.names[subst.idx], names(keys))
  for (i in seq_along(subst.idx)) {
    if (!is.na(matched[i]))
      expr[[subst.idx[i]]] <- keys[[matched[i]]]
  }
  expr
}

.join.exprs <- function(exprs) {
  ret <- quote({})
  for (i in seq_along(exprs)) {
    ret[i + 1] <- exprs[[i]]
  }
  ret
}

