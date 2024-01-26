# The Box-Cox hyperbolic distribution --------------------------------------------------------------

# Density
dbchp <- function(x, mu, sigma, lambda, nu, log = FALSE){

  if(is.matrix(x)) d <- ncol(x) else d <- 1L

  maxl <- max(c(length(x), length(mu), length(sigma), length(lambda), length(nu)))

  x <- rep(x, length.out = maxl)
  mu <- rep(mu, length.out = maxl)
  sigma <- rep(sigma, length.out = maxl)
  lambda <- rep(lambda, length.out = maxl)
  #nu <- rep(nu, length.out = maxl)

  pmf <- rep(-Inf, maxl)

  # NaN index
  pmf[which(mu <= 0 | sigma <= 0 | nu <= 0)] <- NaN

  # Positive density index
  id1 <- which(x > 0 & lambda != 0 & !is.nan(pmf))
  id2 <- which(x > 0 & lambda == 0 & !is.nan(pmf))

  # Extended Box-Cox transformation
  z <- vector()

  z[id1] <- ((x[id1] / mu[id1])^lambda[id1] - 1) / (sigma[id1] * lambda[id1])
  z[id2] <- log(x[id2] / mu[id2]) / sigma[id2]

  # pmf[id1] <- (lambda[id1] - 1) * log(x[id1]) + dhyp(z[id1], nu = nu, log = TRUE) -
  #   phyp(1/(sigma[id1] * abs(lambda[id1])), nu = nu, log.p = TRUE) -
  #   lambda[id1] * log(mu[id1]) - log(sigma[id1])

  # Generating distribution
  log_r <- function(u){
    if (is.vector(u))
      u <- matrix(u, nrow = length(u))

    n <- dim(u)[1]
    d <- dim(u)[2]

    pmf <- matrix(0, n, d)

    # NaN index
    NaNid <- which(nu <= 0, arr.ind = TRUE)
    pmf[NaNid] <- NaN

    # Positive density index
    id <- which(u >= 0 & !is.nan(pmf), arr.ind = TRUE)
    pmf[id] <- log(0.5) -nu * sqrt(1 + u[id]) - log(besselK(x = nu, nu = 1))

    if(d == 1L) as.numeric(pmf) else pmf
  }

  W <- distr::AbscontDistribution(d = function(x) exp(log_r(x^2)))

  pmf[id1] <- (lambda[id1] - 1) * log(x[id1]) + dhyp(z[id1], nu = nu, log = TRUE) -
    log(distr::p(W)(1 / (sigma[id1] * abs(lambda[id1])))) -
    lambda[id1] * log(mu[id1]) - log(sigma[id1])

  pmf[id2] <- dhyp(z[id2], nu = nu, log = TRUE) - log(sigma[id2] * x[id2])

  if(!log) pmf <- exp(pmf)

  if(d > 1L) matrix(pmf, ncol = d) else pmf

}

## Distribution function
pbchp <- function(q, mu, sigma, lambda, nu, lower.tail = TRUE){

  if(is.matrix(q)) d <- ncol(q) else d <- 1L

  maxl <- max(c(length(q), length(mu), length(sigma), length(lambda), length(nu)))

  q <- rep(q, length.out = maxl)
  mu <- rep(mu, length.out = maxl)
  sigma <- rep(sigma, length.out = maxl)
  lambda <- rep(lambda, length.out = maxl)
  #nu <- rep(nu, length.out = maxl)

  # Extended Box-Cox transformation
  z <- rep(NaN, length.out = maxl)

  id1 <- which(q > 0 & mu > 0 & sigma > 0 & lambda != 0 & nu > 0)
  id2 <- which(q > 0 & mu > 0 & sigma > 0 & lambda == 0 & nu > 0)

  z[id1] <- ((q[id1] / mu[id1])^lambda[id1] - 1) / (sigma[id1] * lambda[id1])
  z[id2] <- log(q[id2] / mu[id2]) / sigma[id2]

  id1 <- which(q > 0 & mu > 0 & sigma > 0 & lambda <= 0 & nu > 0)
  id2 <- which(q > 0 & mu > 0 & sigma > 0 & lambda > 0 & nu > 0)

  cdf <- rep(NaN, length.out = maxl)
  # cdf[id1] <- phyp(z[id1], nu = nu[id1]) / phyp(1 / (sigma[id1] * abs(lambda[id1])), nu = nu[id1])
  # cdf[id2] <- (phyp(z[id2], nu = nu[id2]) -  phyp(- 1 / (sigma[id2] * lambda[id2]), nu = nu[id2])) /
  #   phyp(1 / (sigma[id2] * lambda[id2]), nu = nu[id2])

  # Generating distribution
  log_r <- function(u){
    if (is.vector(u))
      u <- matrix(u, nrow = length(u))

    n <- dim(u)[1]
    d <- dim(u)[2]

    pmf <- matrix(0, n, d)

    # NaN index
    NaNid <- which(nu <= 0, arr.ind = TRUE)
    pmf[NaNid] <- NaN

    # Positive density index
    id <- which(u >= 0 & !is.nan(pmf), arr.ind = TRUE)
    pmf[id] <- log(0.5) -nu * sqrt(1 + u[id]) - log(besselK(x = nu, nu = 1))

    if(d == 1L) as.numeric(pmf) else pmf
  }

  W <- distr::AbscontDistribution(d = function(x) exp(log_r(x^2)))

  cdf[id1] <- distr::p(W)(z[id1]) / distr::p(W)(1 / (sigma[id1] * abs(lambda[id1])))
  cdf[id2] <- (distr::p(W)(z[id2]) - distr::p(W)(- 1 / (sigma[id2] * lambda[id2]))) /
    distr::p(W)(1 / (sigma[id2] * lambda[id2]))

  cdf[which(q <= 0 & mu > 0 & sigma > 0 & nu > 0)] <- 0

  if (!lower.tail) cdf <- 1 - cdf

  if(d > 1L) matrix(cdf, ncol = d) else cdf

}

## Quantile function
qbchp <- function(p, mu, sigma, lambda, nu, lower.tail = TRUE){

  if(is.matrix(p)) d <- ncol(p) else d <- 1L

  maxl <- max(c(length(p), length(mu), length(sigma), length(lambda), length(nu)))

  p <- rep(p, length.out = maxl)
  mu <- rep(mu, length.out = maxl)
  sigma <- rep(sigma, length.out = maxl)
  lambda <- rep(lambda, length.out = maxl)
  #nu <- rep(nu, length.out = maxl)

  if (!lower.tail) p <- 1 - p

  qtf <- zp <- rep(NaN, length.out = maxl)

  # z_p
  id1 <- which(p > 0 & p < 1 & mu > 0 & sigma > 0 & lambda <= 0 & nu > 0)
  id2 <- which(p > 0 & p < 1 & mu > 0 & sigma > 0 & lambda > 0 & nu > 0)

  # Generating distribution
  log_r <- function(u){
    if (is.vector(u))
      u <- matrix(u, nrow = length(u))

    n <- dim(u)[1]
    d <- dim(u)[2]

    pmf <- matrix(0, n, d)

    # NaN index
    NaNid <- which(nu <= 0, arr.ind = TRUE)
    pmf[NaNid] <- NaN

    # Positive density index
    id <- which(u >= 0 & !is.nan(pmf), arr.ind = TRUE)
    pmf[id] <- log(0.5) -nu * sqrt(1 + u[id]) - log(besselK(x = nu, nu = 1))

    if(d == 1L) as.numeric(pmf) else pmf
  }

  W <- distr::AbscontDistribution(d = function(x) exp(log_r(x^2)))

  zp[id1] <- distr::q(W)(p[id1] * distr::p(W)(1 / (sigma[id1] * abs(lambda[id1]))))
  zp[id2] <- distr::q(W)(1 - (1 - p[id2]) * distr::p(W)(1 / (sigma[id2] * abs(lambda[id2]))))

  # Quantile function
  id1 <- which(p > 0 & p < 1 & mu > 0 & sigma > 0 & lambda != 0 & nu > 0)
  id2 <- which(p > 0 & p < 1 & mu > 0 & sigma > 0 & lambda == 0 & nu > 0)
  id3 <- which(p == 0 & mu > 0 & sigma > 0 & nu > 0)
  id4 <- which(p == 1 & mu > 0 & sigma > 0 & nu > 0)

  qtf[id1] <- exp(log(mu[id1]) + (1 / lambda[id1]) * log1p(sigma[id1] * lambda[id1] * zp[id1]))
  qtf[id2] <- exp(log(mu[id2]) + sigma[id2] * zp[id2])
  qtf[id3] <- 0
  qtf[id4] <- Inf

  if(d > 1L) matrix(qtf, ncol = d) else qtf
}

# Random generation
rbchp <- function(n, mu, sigma, lambda, nu){
  u <- stats::runif(n)
  qbchp(u, mu, sigma, lambda, nu)
}

# BCS class
bchp <- function(x) {
  out <- list()

  # Abbreviation
  out$abb <- "bchp"

  # Name
  out$name <- "Box-Cox Hyperbolic"

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

# The Log-hyperbolic distribution -------------------------------------------------------------------

# Density
dlhp <- function(x, mu, sigma, nu, log = FALSE, ...) {
  if (is.matrix(x)) d <- ncol(x) else d <- 1L

  maxl <- max(c(length(x), length(mu), length(sigma), length(nu)))

  x <- rep(x, length.out = maxl)
  mu <- rep(mu, length.out = maxl)
  sigma <- rep(sigma, length.out = maxl)

  pmf <- rep(-Inf, maxl)

  # NaN index
  pmf[which(mu <= 0 | sigma <= 0 | nu <= 0)] <- NaN

  # Positive density index
  id <- which(x > 0 & !is.nan(pmf))

  # Transformations
  z <- rep(NaN, length.out = maxl)

  z[id] <- log(x[id] / mu[id]) / sigma[id]


  pmf[id] <- dhyp(z[id], nu = nu, log = TRUE) - log(sigma[id] * x[id])

  if (!log) pmf <- exp(pmf)

  if (d > 1L) matrix(pmf, ncol = d) else pmf
}

## Distribution function
plhp <- function(q, mu, sigma, nu, lower.tail = TRUE, ...) {
  if (is.matrix(q)) d <- ncol(q) else d <- 1L

  maxl <- max(c(length(q), length(mu), length(sigma), length(nu)))

  q <- rep(q, length.out = maxl)
  mu <- rep(mu, length.out = maxl)
  sigma <- rep(sigma, length.out = maxl)

  # Extended Box-Cox transformation
  z <- rep(NaN, length.out = maxl)

  id <- which(q > 0 & mu > 0 & sigma > 0 & nu > 0)

  z[id] <- log(q[id] / mu[id]) / sigma[id]

  cdf <- rep(NaN, length.out = maxl)

  # Generating distribution
  log_r <- function(u){
    if (is.vector(u))
      u <- matrix(u, nrow = length(u))

    n <- dim(u)[1]
    d <- dim(u)[2]

    pmf <- matrix(0, n, d)

    # NaN index
    NaNid <- which(nu <= 0, arr.ind = TRUE)
    pmf[NaNid] <- NaN

    # Positive density index
    id <- which(u >= 0 & !is.nan(pmf), arr.ind = TRUE)
    pmf[id] <- log(0.5) -nu * sqrt(1 + u[id]) - log(besselK(x = nu, nu = 1))

    if(d == 1L) as.numeric(pmf) else pmf
  }

  W <- distr::AbscontDistribution(d = function(x) exp(log_r(x^2)))

  # cdf[id] <- phyp(z[id], nu = nu)
  cdf[id] <- distr::p(W)(z[id])

  cdf[which(q <= 0 & mu > 0 & sigma > 0 & nu > 0)] <- 0

  if (!lower.tail) cdf <- 1 - cdf

  if (d > 1L) matrix(cdf, ncol = d) else cdf
}

## Quantile function
qlhp <- function(p, mu, sigma, nu, lower.tail = TRUE, ...) {
  if (is.matrix(p)) d <- ncol(p) else d <- 1L

  maxl <- max(c(length(p), length(mu), length(sigma), length(nu)))

  p <- rep(p, length.out = maxl)
  mu <- rep(mu, length.out = maxl)
  sigma <- rep(sigma, length.out = maxl)

  if (!lower.tail) p <- 1 - p

  qtf <- zp <- rep(NaN, length.out = maxl)

  # z_p
  id <- which(p > 0 & p < 1 & mu > 0 & sigma > 0 & nu > 0)

  # Generating distribution
  log_r <- function(u){
    if (is.vector(u))
      u <- matrix(u, nrow = length(u))

    n <- dim(u)[1]
    d <- dim(u)[2]

    pmf <- matrix(0, n, d)

    # NaN index
    NaNid <- which(nu <= 0, arr.ind = TRUE)
    pmf[NaNid] <- NaN

    # Positive density index
    id <- which(u >= 0 & !is.nan(pmf), arr.ind = TRUE)
    pmf[id] <- log(0.5) -nu * sqrt(1 + u[id]) - log(besselK(x = nu, nu = 1))

    if(d == 1L) as.numeric(pmf) else pmf
  }

  W <- distr::AbscontDistribution(d = function(x) exp(log_r(x^2)))

  zp[id] <- distr::q(W)(p[id])

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
rlhp <- function(n, mu, sigma, nu) {
  exp(log(mu) + sigma * rhyp(n, nu = nu))
}

# BCS class
lhp <- function(x) {
  out <- list()

  # Abbreviation
  out$abb <- "lhp"

  # Name
  out$name <- "Log-Hyperbolic"

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



# The hyperbolic distribution ----------------------------------------------------------------------

## Probability density function
dhyp <- function(x, mu = 0, sigma = 1, nu, log = FALSE) {
  if (is.matrix(x)) d <- ncol(x) else d <- 1L

  maxl <- max(c(length(x), length(mu), length(sigma), length(nu)))

  x <- rep(x, length.out = maxl)
  mu <- rep(mu, length.out = maxl)
  sigma <- rep(sigma, length.out = maxl)
  nu <- rep(nu, length.out = maxl)

  pmf <- rep(-Inf, length.out = maxl)

  # NaN index
  pmf[which(sigma <= 0 | nu <= 0)] <- NaN

  id <- which(!is.nan(pmf))
  pmf[id] <- hyp_PDF(x[id], mu[id], sigma[id], nu[id], log)

  if (d > 1L) matrix(pmf, ncol = d) else pmf
}

## Cumulative distribution function
phyp <- function(q, mu = 0, sigma = 1, nu, log.p = FALSE) {
  if (is.matrix(q)) d <- ncol(q) else d <- 1L

  maxl <- max(c(length(q), length(mu), length(sigma), length(nu)))

  q <- rep(q, length.out = maxl)
  mu <- rep(mu, length.out = maxl)
  sigma <- rep(sigma, length.out = maxl)

  cdf <- rep(0, length.out = maxl)

  # NaN index
  cdf[which(sigma <= 0 | nu <= 0)] <- NaN

  # Positive density index
  id1 <- which(is.finite(q) & q != mu & !is.nan(cdf))
  id2 <- which(q == mu & !is.nan(cdf))
  id3 <- which(q == -Inf)
  id4 <- which(q == Inf)

  # Constructing the Hyperbolic distribution
  W <- distr::AbscontDistribution(
    d = function(x) dhyp(x, nu = nu),
    Symmetry = distr::SphericalSymmetry(0)
  )

  cdf[id1] <- distr::p(W)((q[id1] - mu[id1]) / sigma[id1])
  cdf[id2] <- 0.5
  cdf[id3] <- 0
  cdf[id4] <- 1


  if (log.p) cdf <- log(cdf)

  if (d > 1L) matrix(cdf, ncol = d) else cdf
}

# Quantile function
qhyp <- function(p, mu = 0, sigma = 1, nu) {
  if (is.matrix(p)) d <- ncol(p) else d <- 1L

  maxl <- max(c(length(p), length(mu), length(sigma), length(nu)))

  p <- rep(p, length.out = maxl)
  mu <- rep(mu, length.out = maxl)
  sigma <- rep(sigma, length.out = maxl)

  qtf <- rep(NA, maxl)

  # NaN index
  qtf[which(p < 0 | p > 1 | sigma <= 0 | nu <= 0)] <- NaN

  # Positive density index
  id1 <- which(p != 0.5 & p > 0 & p < 1 & !is.nan(qtf), arr.ind = TRUE)
  id2 <- which(p == 0 & !is.nan(qtf), arr.ind = TRUE)
  id3 <- which(p == 0.5 & !is.nan(qtf), arr.ind = TRUE)
  id4 <- which(p == 1 & !is.nan(qtf), arr.ind = TRUE)

  # Constructing the Hyperbolic distribution
  W <- distr::AbscontDistribution(
    d = function(x) dhyp(x, nu = nu),
    Symmetry = distr::SphericalSymmetry(0)
  )

  qtf[id1] <- distr::q(W)(p[id1])
  qtf[id2] <- -Inf
  qtf[id3] <- 0
  qtf[id4] <- Inf

  qtf <- mu + sigma * qtf

  if (d > 1L) matrix(qtf, ncol = d) else qtf
}

# Random generation
rhyp <- function(n, mu = 0, sigma = 1, nu) {
  mu + sigma * stats::rnorm(n) * sqrt(GIGrvg::rgig(n, 1.0, 1.0, nu^2))
}

