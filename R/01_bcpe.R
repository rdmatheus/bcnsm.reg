# The Box-Cox power exponential distribution -------------------------------------------------------

# Density
dbcpe <- function(x, mu, sigma, lambda, nu, log = FALSE) {
  if (is.matrix(x)) d <- ncol(x) else d <- 1L

  maxl <- max(c(length(x), length(mu), length(sigma), length(lambda), length(nu)))

  x <- rep(x, length.out = maxl)
  mu <- rep(mu, length.out = maxl)
  sigma <- rep(sigma, length.out = maxl)
  lambda <- rep(lambda, length.out = maxl)
  nu <- rep(nu, length.out = maxl)

  pmf <- rep(-Inf, maxl)

  # NaN index
  pmf[which(mu <= 0 | sigma <= 0 | nu < 1)] <- NaN

  # Positive density index
  id1 <- which(x > 0 & lambda != 0 & !is.nan(pmf))
  id2 <- which(x > 0 & lambda == 0 & !is.nan(pmf))

  # Extended Box-Cox transformation
  z <- rep(NaN, length.out = maxl)

  z[id1] <- ((x[id1] / mu[id1])^lambda[id1] - 1) / (sigma[id1] * lambda[id1])
  z[id2] <- log(x[id2] / mu[id2]) / sigma[id2]


  pmf[id1] <- (lambda[id1] - 1) * log(x[id1]) + rPE(z[id1]^2, nu = nu[id1], log = TRUE) -
    log(RPE(1 / (sigma[id1] * abs(lambda[id1])), nu = nu[id1])) -
    lambda[id1] * log(mu[id1]) - log(sigma[id1])

  pmf[id2] <- rPE(z[id2]^2, nu = nu[id2], log = TRUE) - log(sigma[id2] * x[id2])

  if (!log) pmf <- exp(pmf)

  if (d > 1L) matrix(pmf, ncol = d) else pmf
}

## Distribution function
pbcpe <- function(q, mu, sigma, lambda, nu, lower.tail = TRUE) {
  if (is.matrix(q)) d <- ncol(q) else d <- 1L

  maxl <- max(c(length(q), length(mu), length(sigma), length(lambda), length(nu)))

  q <- rep(q, length.out = maxl)
  mu <- rep(mu, length.out = maxl)
  sigma <- rep(sigma, length.out = maxl)
  lambda <- rep(lambda, length.out = maxl)
  nu <- rep(nu, length.out = maxl)

  # Extended Box-Cox transformation
  z <- rep(NaN, length.out = maxl)

  id1 <- which(q > 0 & mu > 0 & sigma > 0 & lambda != 0 & nu >= 1)
  id2 <- which(q > 0 & mu > 0 & sigma > 0 & lambda == 0 & nu >= 1)

  z[id1] <- ((q[id1] / mu[id1])^lambda[id1] - 1) / (sigma[id1] * lambda[id1])
  z[id2] <- log(q[id2] / mu[id2]) / sigma[id2]

  id1 <- which(q > 0 & mu > 0 & sigma > 0 & lambda <= 0 & nu >= 1)
  id2 <- which(q > 0 & mu > 0 & sigma > 0 & lambda > 0 & nu >= 1)

  cdf <- rep(NaN, length.out = maxl)
  cdf[id1] <- RPE(z[id1], nu = nu[id1]) / RPE(1 / (sigma[id1] * abs(lambda[id1])), nu = nu[id1])
  cdf[id2] <- (RPE(z[id2], nu = nu[id2]) - RPE(-1 / (sigma[id2] * lambda[id2]), nu = nu[id2])) /
    RPE(1 / (sigma[id2] * lambda[id2]), nu = nu[id2])

  cdf[which(q <= 0 & mu > 0 & sigma > 0 & nu >= 1)] <- 0

  if (!lower.tail) cdf <- 1 - cdf

  if (d > 1L) matrix(cdf, ncol = d) else cdf
}

## Quantile function
qbcpe <- function(p, mu, sigma, lambda, nu, lower.tail = TRUE) {
  if (is.matrix(p)) d <- ncol(p) else d <- 1L

  maxl <- max(c(length(p), length(mu), length(sigma), length(lambda), length(nu)))

  p <- rep(p, length.out = maxl)
  mu <- rep(mu, length.out = maxl)
  sigma <- rep(sigma, length.out = maxl)
  lambda <- rep(lambda, length.out = maxl)
  nu <- rep(nu, length.out = maxl)

  if (!lower.tail) p <- 1 - p

  qtf <- zp <- rep(NaN, length.out = maxl)

  # z_p
  id1 <- which(p > 0 & p < 1 & mu > 0 & sigma > 0 & lambda <= 0 & nu >= 1)
  id2 <- which(p > 0 & p < 1 & mu > 0 & sigma > 0 & lambda > 0 & nu >= 1)

  zp[id1] <- qPE(p[id1] * RPE(1 / (sigma[id1] * abs(lambda[id1])), nu = nu[id1]), nu = nu[id1])
  zp[id2] <- qPE(1 - (1 - p[id2]) * RPE(1 / (sigma[id2] * abs(lambda[id2])), nu = nu[id2]), nu = nu[id2])

  # Quantile function
  id1 <- which(p > 0 & p < 1 & mu > 0 & sigma > 0 & lambda != 0 & nu >= 1)
  id2 <- which(p > 0 & p < 1 & mu > 0 & sigma > 0 & lambda == 0 & nu >= 1)
  id3 <- which(p == 0 & mu > 0 & sigma > 0 & nu >= 1)
  id4 <- which(p == 1 & mu > 0 & sigma > 0 & nu >= 1)

  qtf[id1] <- exp(log(mu[id1]) + (1 / lambda[id1]) * log1p(sigma[id1] * lambda[id1] * zp[id1]))
  qtf[id2] <- exp(log(mu[id2]) + sigma[id2] * zp[id2])
  qtf[id3] <- 0
  qtf[id4] <- Inf

  if (d > 1L) matrix(qtf, ncol = d) else qtf
}

# Random generation
rbcpe <- function(n, mu, sigma, lambda, nu) {
  u <- stats::runif(n)
  qbcpe(u, mu, sigma, lambda, nu)
}

# BCS class
bcpe <- function(x) {
  out <- list()

  # Abbreviation
  out$abb <- "bcpe"

  # Name
  out$name <- "Box-Cox Power Exponential"

  # Number of parameters
  out$npar <- 4

  # Extra parameter
  out$extrap <- TRUE

  # Initial values -------------------------------------------------------------
  out$start <- function(x) {

    n <- length(x)

    gamlss_fit <- suppressWarnings(try(gamlss::gamlss(x ~ 1, family = gamlss.dist::BCPE(mu.link = "log"), trace = FALSE), silent = TRUE))

    if (unique(grepl("Error", gamlss_fit))) {
      convergence <- FALSE
    } else {
      convergence <- gamlss_fit$converged
    }

    if (convergence) {
      c(exp(stats::coef(gamlss_fit, "mu")), exp(stats::coef(gamlss_fit, "sigma")),
        stats::coef(gamlss_fit, "nu"), min(exp(stats::coef(gamlss_fit, "tau")), 20))
    } else {

      CV <- 0.75 * (stats::quantile(x, 0.75) - stats::quantile(x, 0.25)) / stats::median(x)

      mu0 <- stats::median(x)
      sigma0 <- asinh(CV / 1.5) * stats::qnorm(0.75)

      z <- log(x / mu0) / sigma0

      grid <- seq(1, 20, 1)
      upsilon <- function(nu){
        cdf <- sort(gamlss.dist::dPE(z, mu = 0, sigma = 1, nu = nu))
        temp <- stats::qqnorm(stats::qnorm(cdf), plot.it = FALSE)
        mean(abs(sort(temp$x) - sort(temp$y)))
      }

      out <- apply(matrix(grid), 1, upsilon)
      nu0 <- grid[which.min(out)]

      c(mu0, sigma0, 0L, nu0)
    }

  }

  structure(out, class = "bcs")
}


# The Log-power exponential distribution -------------------------------------------------------------------

# Density
dlpe <- function(x, mu, sigma, nu, log = FALSE, ...) {
  if (is.matrix(x)) d <- ncol(x) else d <- 1L

  maxl <- max(c(length(x), length(mu), length(sigma), length(nu)))

  x <- rep(x, length.out = maxl)
  mu <- rep(mu, length.out = maxl)
  sigma <- rep(sigma, length.out = maxl)
  nu <- rep(nu, length.out = maxl)

  pmf <- rep(-Inf, maxl)

  # NaN index
  pmf[which(mu <= 0 | sigma <= 0 | nu <= 0)] <- NaN

  # Positive density index
  id <- which(x > 0 & !is.nan(pmf))

  # Transformations
  z <- rep(NaN, length.out = maxl)

  z[id] <- log(x[id] / mu[id]) / sigma[id]

  pmf[id] <- rPE(z[id]^2, nu = nu[id], log = TRUE) - log(sigma[id] * x[id])

  if (!log) pmf <- exp(pmf)

  if (d > 1L) matrix(pmf, ncol = d) else pmf
}

## Distribution function
plpe <- function(q, mu, sigma, nu, lower.tail = TRUE, ...) {
  if (is.matrix(q)) d <- ncol(q) else d <- 1L

  maxl <- max(c(length(q), length(mu), length(sigma), length(nu)))

  q <- rep(q, length.out = maxl)
  mu <- rep(mu, length.out = maxl)
  sigma <- rep(sigma, length.out = maxl)
  nu <- rep(nu, length.out = maxl)

  # Extended Box-Cox transformation
  z <- rep(NaN, length.out = maxl)

  id <- which(q > 0 & mu > 0 & sigma > 0 & nu > 0)

  z[id] <- log(q[id] / mu[id]) / sigma[id]

  cdf <- rep(NaN, length.out = maxl)
  cdf[id] <- RPE(z[id], nu = nu[id])

  cdf[which(q <= 0 & mu > 0 & sigma > 0 & nu > 0)] <- 0

  if (!lower.tail) cdf <- 1 - cdf

  if (d > 1L) matrix(cdf, ncol = d) else cdf
}

## Quantile function
qlpe <- function(p, mu, sigma, nu, lower.tail = TRUE, ...) {
  if (is.matrix(p)) d <- ncol(p) else d <- 1L

  maxl <- max(c(length(p), length(mu), length(sigma), length(nu)))

  p <- rep(p, length.out = maxl)
  mu <- rep(mu, length.out = maxl)
  sigma <- rep(sigma, length.out = maxl)
  nu <- rep(nu, length.out = maxl)

  if (!lower.tail) p <- 1 - p

  qtf <- zp <- rep(NaN, length.out = maxl)

  # z_p
  id <- which(p > 0 & p < 1 & mu > 0 & sigma > 0 & nu > 0)

  zp[id] <- qPE(p[id], nu = nu[id])

  # Quantile function
  id1 <- which(p > 0 & p < 1 & mu > 0 & sigma > 0 & nu > 0)
  id2 <- which(p == 0 & mu > 0 & sigma > 0 & nu > 0)
  id3 <- which(p == 1 & mu > 0 & sigma > 0 & nu > 0)

  qtf[id1] <- exp(log(mu[id1]) + sigma[id1] * zp[id1])
  qtf[id2] <- 0
  qtf[id3] <- Inf

  if (d > 1L) matrix(qtf, ncol = d) else qtf
}

# Random generation
rlpe <- function(n, mu, sigma, nu) {
  u <- stats::runif(n)
  exp(log(mu) + sigma * qPE(u, nu = nu))
}

# BCS class
lpe <- function(x) {
  out <- list()

  # Abbreviation
  out$abb <- "lpe"

  # Name
  out$name <- "Log-Power Exponential"

  # Number of parameters
  out$npar <- 3

  # Extra parameter
  out$extrap <- TRUE

  # Initial values -------------------------------------------------------------
  out$start <- function(x) {

    n <- length(x)

    gamlss_fit <- suppressWarnings(try(gamlss::gamlss(x ~ 1, family = gamlss.dist::BCPE(mu.link = "log"),
                                                      trace = FALSE, nu.fix = TRUE, nu.start = 0L), silent = TRUE))

    if (unique(grepl("Error", gamlss_fit))) {
      convergence <- FALSE
    } else {
      convergence <- gamlss_fit$converged
    }

    if (convergence) {
      c(exp(stats::coef(gamlss_fit, "mu")), exp(stats::coef(gamlss_fit, "sigma")),
        min(exp(stats::coef(gamlss_fit, "tau")), 20))
    } else {
      CV <- 0.75 * (stats::quantile(x, 0.75) - stats::quantile(x, 0.25)) / stats::median(x)

      mu0 <- stats::median(x)
      sigma0 <- asinh(CV / 1.5) * stats::qnorm(0.75)

      z <- log(x / mu0) / sigma0

      grid <- seq(1, 20, 1)
      upsilon <- function(nu){
        cdf <- sort(gamlss.dist::dPE(z, mu = 0, sigma = 1, nu = nu))
        temp <- stats::qqnorm(stats::qnorm(cdf), plot.it = FALSE)
        mean(abs(sort(temp$x) - sort(temp$y)))
      }

      out <- apply(matrix(grid), 1, upsilon)
      nu0 <- grid[which.min(out)]

      c(mu0, sigma0, nu0)

    }

  }


  structure(out, class = "bcs")
}

# The standard power exponential distribution ------------------------------------------------------

## Density generating function
rPE <- function(u, nu, log = FALSE) {
  maxl <- max(c(length(u), length(nu)))

  u <- rep(u, length.out = maxl)
  nu <- rep(nu, length.out = maxl)

  pmf <- rep(-Inf, length.out = maxl)

  # NaN index
  pmf[which(nu < 1)] <- NaN

  # Positive density index
  id <- which(u >= 0L & !is.nan(pmf))

  log.c <- 0.5 * (-(2 / nu[id]) * log(2) + lgamma(1 / nu[id]) - lgamma(3 / nu[id]))
  c <- exp(log.c)
  pmf[id] <- log(nu[id]) - log.c - (1 + 1 / nu[id]) * log(2) - lgamma(1 / nu[id]) -
    0.5 * (u[id]^(nu[id] / 2)) / (c^nu[id])

  if (!log) pmf <- exp(pmf)

  pmf
}

## Cumulative distribution function
RPE <- function(q, nu) {
  maxl <- max(c(length(q), length(nu)))

  q <- rep(q, length.out = maxl)
  nu <- rep(nu, length.out = maxl)

  cdf <- rep(0, length.out = maxl)

  # NaN index
  cdf[which(nu < 1)] <- NaN

  # Positive density index
  id1 <- which(is.finite(q) & q != 0 & nu < 10000 & !is.nan(cdf))
  id1_aux <- which(is.finite(q) & q != 0 & nu >= 10000 & !is.nan(cdf))
  id2 <- which(q == 0 & !is.nan(cdf))
  id3 <- which(q == -Inf)
  id4 <- which(q == Inf)

  log.c <- 0.5 * (-(2 / nu[id1]) * log(2) + lgamma(1 / nu[id1]) - lgamma(3 / nu[id1]))
  c <- exp(log.c)
  s <- 0.5 * ((abs(q[id1] / c))^nu[id1])

  cdf[id1] <- 0.5 * (1 + stats::pgamma(s, shape = 1 / nu[id1], scale = 1) * sign(q[id1]))
  cdf[id1_aux] <- (q[id1_aux] + sqrt(3)) / sqrt(12)
  cdf[id2] <- 0.5
  cdf[id3] <- 0
  cdf[id4] <- 1

  cdf
}

# Quantile function
qPE <- function(p, nu) {
  maxl <- max(c(length(p), length(nu)))

  p <- rep(p, length.out = maxl)
  nu <- rep(nu, length.out = maxl)

  qtf <- rep(NA, maxl)

  # NaN index
  qtf[which(nu < 1)] <- NaN

  # Positive density index
  id1 <- which(p != 0.5 & p > 0 & p < 1 & !is.nan(qtf), arr.ind = TRUE)
  id2 <- which(p == 0 & !is.nan(qtf), arr.ind = TRUE)
  id3 <- which(p == 0.5 & !is.nan(qtf), arr.ind = TRUE)
  id4 <- which(p == 1 & !is.nan(qtf), arr.ind = TRUE)

  log.c <- 0.5 * (-(2 / nu[id1]) * log(2) + lgamma(1 / nu[id1]) - lgamma(3 / nu[id1]))
  c <- exp(log.c)
  s <- stats::qgamma((2 * p[id1] - 1) * sign(p[id1] - 0.5), shape = (1 / nu[id1]), scale = 1)

  qtf[id1] <- sign(p[id1] - 0.5) * ((2 * s)^(1 / nu[id1])) * c
  qtf[id2] <- -Inf
  qtf[id3] <- 0
  qtf[id4] <- Inf

  qtf
}
