# studi con diverso design
#                    lista [2]          lista [2]
#                      --v--               --v--
alpha.cf1 <-  function(beta, mu0, sigma20, sigma2, design) {
  gamma <- get.matrix.from.design(design, "gamma")
  t <- get.matrix.from.design(design, "theta")
  return(c(
    alpha1 = -(((t[1, 1] * beta[1] - t[1, 2]) * sigma20 +
                  (mu0 * beta[1] - t[1, 2]) * gamma[1, 1]) /
                 (sigma20 + gamma[1, 1])),
    alpha2 = -(((t[2, 1] * beta[2] - t[2, 3]) * sigma20 +
                  mu0 * beta[2] * gamma[2, 1] - t[2, 3] * gamma[2, 1]) /
                 (sigma20 + gamma[2, 1]))
  ))
}

sigma.cf1 <- function(beta, mu0, sigma20, design) {
  gamma <- get.matrix.from.design(design, "gamma")
  t <- get.matrix.from.design(design, "theta")
  return(c(
    sigma21 = -(((gamma[1, 2] + beta[1]^2 * gamma[1, 1]) * sigma20 +
                   gamma[1, 1] * gamma[1, 2]) /
                  (sigma20 + gamma[1, 1])),
    sigma22 = -(((gamma[2, 3] + beta[2]^2 * gamma[2, 1]) * sigma20 +
                   gamma[2, 1] * gamma[2, 3]) /
                  (sigma20 + gamma[2, 1]))
  ))
}

# studi con lo stesso design

alpha.cf2 <-  function(beta, mu0, sigma20, sigma2, design) {
  gamma <- get.matrix.from.design(design, "gamma")
  t <- get.matrix.from.design(design, "theta")
  return(
    -((((t[1, 1] * beta[1] - t[1, 2]) * sigma20^2 +
          ((t[1, 1] * beta[1] - t[1, 2]) * gamma[2, 1] +
             (mu0 * beta[1] - t[1, 2]) * gamma[1, 1]) * sigma20 +
          (mu0 * beta[1] - t[1, 2]) * gamma[1, 1] * gamma[2, 1]) * sigma2[1] +
         (t[2, 1] * beta[1] - t[2, 2]) * sigma20^3 +
         ((t[1, 1] * beta[1] - t[1, 2]) * gamma[2, 2] +
            (t[1, 1] * beta[1]^3 - t[1, 2] * beta[1]^2 +
               mu0 * beta[1] - t[2, 2]) * gamma[2, 1] +
            (t[2, 1] * beta[1] - t[2, 2]) * gamma[1, 2] +
            (t[2, 1] * beta[1]^3 - t[2, 2] * beta[1]^2 +
               t[2, 1] * beta[1] - t[2, 2]) * gamma[1, 1]) * sigma20^2 +
         (((t[1, 1] * beta[1] - t[1, 2]) * gamma[2, 1] +
             (mu0 * beta[1] - t[1, 2]) * gamma[1, 1]) * gamma[2, 2] +
            ((mu0 * beta[1] - t[2, 2]) * gamma[1, 2] +
               (2 * mu0 * beta[1]^3 + ( - t[2, 2] - t[1, 2]) * beta[1]^2 +
                  mu0 * beta[1] - t[2, 2]) * gamma[1, 1]) * gamma[2, 1] +
            (t[2, 1] * beta[1] - t[2, 2]) * gamma[1, 1] * gamma[1, 2]) * sigma20 +
         (mu0 * beta[1] - t[1, 2]) * gamma[1, 1] * gamma[2, 1] * gamma[2, 2] +
         (mu0 * beta[1] - t[2, 2]) * gamma[1, 1] * gamma[1, 2] * gamma[2, 1]) /
        ((sigma20^2 + (gamma[2, 1] + gamma[1, 1]) * sigma20 +
            gamma[1, 1] * gamma[2, 1]) * sigma2[1] + sigma20^3 +
           (gamma[2, 2] + (beta[1]^2 + 1) * gamma[2, 1] + gamma[1, 2] +
              (beta[1]^2 + 1) * gamma[1, 1]) * sigma20^2 +
           ((gamma[2, 1] + gamma[1, 1]) * gamma[2, 2] +
              (gamma[1, 2] + (2 * beta[1]^2 + 1) * gamma[1, 1]) * gamma[2, 1] +
              gamma[1, 1] * gamma[1, 2]) * sigma20 +
           gamma[1, 1] * gamma[2, 1] * gamma[2, 2] +
           gamma[1, 1] * gamma[1, 2] * gamma[2, 1]))
  )
}

sigma.cf2 <- function(alpha, beta, mu0, sigma20, design) {
  gamma <- get.matrix.from.design(design, "gamma")
  t <- get.matrix.from.design(design, "theta")
  return(
    -(((gamma[2, 2] + beta[1]^2 * gamma[2, 1] - t[2, 1]^2 * beta[1]^2 +
          (2 * t[2, 1] * t[2, 2] - 2 * t[2, 1] * alpha[1]) * beta[1] - alpha[1]^2 +
          2 * t[2, 2] * alpha[1] - t[2, 2]^2) * sigma20^2 +
         (2 * gamma[2, 1] * gamma[2, 2] + beta[1]^2 * gamma[2, 1]^2 +
            ( - (2 * t[2, 1] * mu0 * beta[1]^2) +
                (( - (2 * mu0) - 2 * t[2, 1]) * alpha[1] + 2 * t[2, 2] * mu0 +
                   2 * t[2, 1] * t[2, 2]) * beta[1] - 2 * alpha[1]^2 +
                  4 * t[2, 2] * alpha[1] - 2 * t[2, 2]^2) * gamma[2, 1]) * sigma20 +
         gamma[2, 1]^2 * gamma[2, 2] + ( -(mu0^2 * beta[1]^2) +
                                           (2 * t[2, 2] * mu0 -
                                              2 * mu0 * alpha[1]) * beta[1] -
                                           alpha[1]^2 +
                                            2 * t[2, 2] * alpha[1] -
                                            t[2, 2]^2) * gamma[2, 1]^2) /
        (sigma20^2 + 2 * gamma[2, 1] * sigma20 + gamma[2, 1]^2))

  )
}
