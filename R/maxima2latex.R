maxima2tex <- function(filenm) {
  sub.list <- c("µ" = "\\mu", "β" = "\\beta", "Σ" = "\\Sigma", "α" = "\\alpha",
                "γ" = "\\gamma",
                "θ" = "\\theta", "ψ" = "\\psi", "Ψ" = "\\Psi", "λ" = "\\lambda",
                "σ" = "\\sigma", "*" = " ", " . " = " ", "log" = "\\log", "exp" = "\\exp")


  txt <- readLines(filenm)
  for (i in seq_along(sub.list)) {
    txt <- gsub(names(sub.list)[i], sub.list[i], txt, fixed = TRUE)
  }

  # pedici e apici
  txt <- stringr::str_replace_all(txt, r"{_([0-9]+)}", r"{_{\1}}")
  txt <- stringr::str_replace_all(txt, r"{\^+\(([^()]*)\)}", r"{^{\1}}")

  # frazioni con parentesi
  txt <- stringr::str_replace_all(txt, r"{\(([^()]*)\) ?/ ?\(([^()]*)\)}", r"{\frac{\1}{\2}}")


  writeLines(txt, stringr::str_replace_all(filenm, "([.]?[a-zA-Z]+?)$", "_tex\\1"))

}
