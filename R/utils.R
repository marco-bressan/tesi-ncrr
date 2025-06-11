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

...not.used <- function(..., output.type = "message") {
  if (length(output.type) == 0)
    invisible()
  if (...length() > 0)
    do.call(output.type, list("Argomenti extra (`...`) non utilizzati: ",
                              paste(...names(), collapse = ", "),
                              "."))
  invisible()
}

##' Effettua una sostituzione nella lista dei parametri di una NCRR
##'
##'
##' @title Sostituisci parametri
##' @param params Lista o vettore dei parametri
##' @param subst Lista o vettore delle sostituzioni.
##' La sostituzione avviene in base al nome, che dunque è obbligatorio per
##' ogni parametro
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
