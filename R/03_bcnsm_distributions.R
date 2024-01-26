dbcnsm <- function(x, mu, sigma, lambda, nu = NULL, Gamma = diag(ncol(x)),
                   copula = c("gaussian", "t", "slash", "hyp"),
                   delta = NULL, margins = "bcno", log = FALSE) {

  ### Setting dimensions
  if (is.vector(x))
    x <- matrix(x, ncol = length(x))

  n <- dim(x)[1]
  d <- dim(x)[2]

  ### Marginal parameters
  mu <- matrix(mu, ncol = d)
  sigma <- matrix(sigma, ncol = d)
  lambda <- matrix(lambda, ncol = d)
  aux_nu <- as.numeric(stats::na.exclude(nu))

  mu <- do.call(rbind, replicate(n / dim(mu)[1], mu, simplify = FALSE))
  sigma <- do.call(rbind, replicate(n / dim(sigma)[1], sigma, simplify = FALSE))
  lambda <- do.call(rbind, replicate(n / dim(lambda)[1], lambda, simplify = FALSE))

  margins <- as.vector(matrix(margins, 1, d))

  nu <- rep(NA, d)
  nu_id <- apply(matrix(margins, ncol = 1), 1, function(x) as.bcs(x)$extrap)
  nu[nu_id] <- aux_nu

  ### Copula
  copula <- match.arg(copula, c("gaussian", "t", "slash", "hyp"))
  copula_inf <- make_copula(copula, delta)

  dgf <- copula_inf$dgf
  qPSI <- copula_inf$qPSI

  EPS <- .Machine$double.eps^(1 / 2)
  w <- log_f <- matrix(NA, n, d)
  for (j in 1:d) {

    w[, j] <- pmin(pmax(get(paste0("p", margins[j]))(q = x[, j],
                                                     mu = mu[, j],
                                                     sigma = sigma[, j],
                                                     lambda = lambda[, j],
                                                     nu = nu[j]), EPS), 1 - EPS)

    log_f[, j] <- get(paste0("d", margins[j]))(x = x[, j],
                                               mu = mu[, j],
                                               sigma = sigma[, j],
                                               lambda = lambda[, j],
                                               nu = nu[j], log = TRUE)

  }

  epsilon <- matrix(qPSI(w), ncol = d)
  dec <- Rfast::cholesky(Gamma)

  tmp <- backsolve(dec, t(epsilon), transpose = TRUE)
  rss <- colSums(tmp^2)

  den <- -sum(log(diag(dec))) + dgf(rss, d, log = TRUE) -
    rowSums(matrix(dgf(epsilon^2, 1L, log = TRUE), ncol = d)) + rowSums(log_f)

  if (log) {
    den
  } else {
    exp(den)
  }
}

# Copula information -------------------------------------------------------------------------------
make_copula <- function(copula, delta) {

  switch(copula,

         gaussian = {

           dPSI <- function(x, log = FALSE) stats::dnorm(x, log = log)
           pPSI <- function(q) stats::pnorm(q)
           qPSI <- function(p) stats::qnorm(p)
           dgf <- function(u, d, log = FALSE) {
             out <- -0.5 * u - (d/2) * log(2 * pi)

             if (!log) out <- exp(out)

             out
           }

         },

         t = {

           dPSI <- function(x, log = FALSE) stats::dt(x, delta, log = log)
           pPSI <- function(q) stats::pt(q, delta)
           qPSI <- function(p) stats::qt(p, delta)
           dgf <- function(u, d, log = FALSE){

             out <- lgamma(0.5 * (delta + d)) - 0.5 * (delta + d) * log(1 + u / delta) -
               lgamma(delta / 2) - (d/2) * log(delta * pi)

             if (!log) out <- exp(out)

             out
           }

         },

         slash = {

           Wsl <- distr::AbscontDistribution(
             d = function(x) dslash(x, nu = delta),
             Symmetry = distr::SphericalSymmetry(0)
           )

           dPSI <- function(x, log = FALSE) dslash(x, nu = delta, log = log)
           pPSI <- function(q) distr::p(Wsl)(q)
           qPSI <- function(p) distr::q(Wsl)(p)
           dgf <- function(u, d, log = FALSE){

             id1 <- which(u == 0)
             id2 <- which(u != 0)

             out <- vector("numeric", length(u))

             out[id1] <- log(2 * delta) - (d/2) * log(2 * pi) - log(2 * delta + d)
             out[id2] <- log(delta) + delta * log(2) + log(ig(delta + d / 2, u[id2] / 2)) -
               (d/2) * log(pi) - (delta + d/2) * log(u[id2])

             if (!log) out <- exp(out)

             out
           }

         },

         hyp = {

           Whp <- distr::AbscontDistribution(
             d = function(x) dhyp(x, nu = delta),
             Symmetry = distr::SphericalSymmetry(0)
           )

           dPSI <- function(x, log = FALSE) dhyp(x, nu = delta, log = log)
           pPSI <- function(q) distr::p(Whp)(q)
           qPSI <- function(p) distr::q(Whp)(p)
           dgf <- function(u, d, log = FALSE){

             out <- (d - 1) * log(delta) + log(besselK(delta * sqrt(1 + u), nu = 1 - d/2)) -
               (d/2) * log(2 * pi) - log(besselK(delta, nu = 1)) -
               (d/2 - 1) * log(delta * sqrt(1 + u))

             if (!log) out <- exp(out)

             out
           }

         },

         stop(gettextf("%s copula not recognised", sQuote(copula)), domain = NA))

  list(dPSI = dPSI, pPSI = pPSI, qPSI = qPSI, dgf = dgf)

}

