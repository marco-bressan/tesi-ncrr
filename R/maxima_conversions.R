maxima2tex <- function(filenm) {
  sub.list <- c("(d/dβ_01)(det(Σ_1)^(-1))" = r"{\left[\frac{\partial}{\partial\beta_{01}}\left|\tilde{\Sigma}_{1}\right|^{-1}\right]}",
                "µ" = "\\mu", "β" = "\\beta", "Σ" = "\\Sigma", "α" = "\\alpha",
                "γ" = "\\gamma",
                "θ" = "\\theta", "ψ" = "\\psi", "Ψ" = "\\Psi", "λ" = "\\lambda",
                "σ" = "\\sigma", "·" = " ", "*" = " ", " . " = " ", "log" = "\\log", "exp" = "\\exp")


  txt <- readLines(filenm)
  for (i in seq_along(sub.list)) {
    txt <- gsub(names(sub.list)[i], sub.list[i], txt, fixed = TRUE)
  }

  # pedici e apici
  txt <- stringr::str_replace_all(txt, r"{_([0-9]+)}", r"{_{\1}}")
  txt <- stringr::str_replace_all(txt, r"{\^+\(([^()]*)\)}", r"{^{\1}}")

  # frazioni con parentesi
  #txt <- stringr::str_replace_all(txt, r"{\(([^()]*)\) ?/ ?\(([^()]*)\)}", r"{\frac{\1}{\2}}")


  writeLines(txt, stringr::str_replace_all(filenm, "([.]?[a-zA-Z]+?)$", "_tex\\1"))

}

maxima2R <- function(filenm) {
  sub.list <- c("µ_0" = "mu0", "σ_0" = "sigma20",
                "µ" = "mu", "β" = "beta", "Σ" = "Sigma", "α" = "alpha",
                "γ" = "gamma",
                "θ" = "theta", "ψ" = "psi", "Ψ" = "Psi", "λ" = "lambda",
                "σ" = "sigma2", "*" = " * ", "·" = " * ",
                "+" = " + ", "-" = " - ", "/" = " / ","." = "%*%")

  txt <- readLines(filenm)
  for (i in seq_along(sub.list)) {
    txt <- gsub(names(sub.list)[i], sub.list[i], txt, fixed = TRUE)
  }
  txt <- stringr::str_replace_all(txt, "^([a-zA-Z '.,;]{10,})", "# \\1")

  txt <- stringr::str_replace_all(txt, "([a-z]+)_0([0-9])", "\\1[\\2]")
  txt <- stringr::str_replace_all(txt, "([a-z]+)_([1-9])([0-9])", "\\1[\\2, \\3]")

  writeLines(txt, gsub("([.]?[a-zA-Z]+?)$", "__.R", filenm))

}
