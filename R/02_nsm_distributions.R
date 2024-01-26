# The Multivariate Normal Scale Mixture Distributions ----------------------------------------------

## Multivariate normal -----------------------------------------------------------------------------

# Probability density function
dmv_gaussian <- function(x, mu, Sigma, delta = NULL, log = FALSE) {

  mvnfast::dmvn(x, mu, Sigma, log = log, ncores = 2)

}

# Random generation
rmv_gaussian <- function(n, mu, Sigma, delta = NULL) {

  d <- ncol(Sigma)
  A <- Rfast::cholesky(Sigma)
  Z <- matrix(stats::rnorm(n * d), ncol = d)

  mu + Z%*%A

}

# Distribution function of the squared Mahalanobis distance
maha_gaussian <- function(q, delta = NULL, d) {
  stats::pchisq(q, d)
}

# Conditional distribution
dcmv_gaussian <- function(x, mu, Sigma, given_id, delta = NULL, log = FALSE) {

  d <- ncol(Sigma)
  d1 <- length(given_id)
  d2 <- d - d1

  mu1 <- mu[given_id]
  mu2 <- mu[-given_id]

  Sigma11 <- matrix(Sigma[1:d1, 1:d1], d1, d1)
  Sigma12 <- matrix(Sigma[1:d1, (d1 + 1):d], d1, d2)
  Sigma21 <- t(Sigma12)
  Sigma22 <- matrix(Sigma[(d1 + 1):d, (d1 + 1):d], d2, d2)

  aux <- function(x) {

    X_given <- x[given_id]

    mu_star <- mu2 + Sigma21%*%solve(Sigma11)%*%(X_given - mu1)
    Sigma_star <- Sigma22 - Sigma21%*%solve(Sigma11)%*%Sigma12

    mvnfast::dmvn(x[-given_id], mu = mu_star, sigma = Sigma_star, log = log, ncores = 2)


  }

  if (is.vector(x)) {
    aux(x)
  } else {
    apply(x, 1, aux)
  }

}


## Multivariate Student's t ------------------------------------------------------------------------

# X = Z/sqrt(W),     Z ~ Nd(mu, Sigma) and W ~ Gamma(delta/2, delta/2)

# Probability density function
dmv_t <- function(x, mu, Sigma, delta = NULL, log = FALSE) {

  mvnfast::dmvt(x, mu, Sigma, log = log, df = delta, ncores = 2)

}

# Random generation
rmv_t <- function(n, mu, Sigma, delta) {

  d <- ncol(Sigma)
  A <- Rfast::cholesky(Sigma)
  Z <- matrix(stats::rnorm(n * d), ncol = d)

  mu + Z%*%A / sqrt(stats::rgamma(n, delta / 2, delta / 2))

}

# Distribution function of the squared Mahalanobis distance
maha_t <- function(q, delta, d) {
  stats::pf(q / d, d, delta)
}

# Conditional distribution
## https://en.wikipedia.org/wiki/Multivariate_t-distribution
dcmv_t <- function(x, mu, Sigma, given_id, delta, log = FALSE) {

  d <- ncol(Sigma)
  d1 <- length(given_id)
  d2 <- d - d1

  mu1 <- mu[given_id]
  mu2 <- mu[-given_id]

  Sigma11 <- matrix(Sigma[1:d1, 1:d1], d1, d1)
  Sigma12 <- matrix(Sigma[1:d1, (d1 + 1):d], d1, d2)
  Sigma21 <- t(Sigma12)
  Sigma22 <- matrix(Sigma[(d1 + 1):d, (d1 + 1):d], d2, d2)


  aux <- function(x) {

    X_given <- x[given_id]

    M <- Rfast::mahala(X_given, mu1, Sigma11)

    mu_star <- mu2 + Sigma21%*%solve(Sigma11)%*%(X_given - mu1)
    Sigma_star <- Sigma22 - Sigma21%*%solve(Sigma11)%*%Sigma12

    mvnfast::dmvt(x[-given_id], mu = mu_star,
                  sigma = ((delta + M) / (delta + d1)) * Sigma_star,
                  df = delta + d1,
                  log = log, ncores = 2)


  }

  if (is.vector(x)) {
    aux(x)
  } else {
    apply(x, 1, aux)
  }

}


## Multivariate slash ------------------------------------------------------------------------------

# X = Z/sqrt(W),     Z ~ Nd(mu, Sigma) and W ~ Beta(delta, 1)

# Probability density function
dmv_slash <- function(x, mu, Sigma, delta, log = FALSE) {

  delta <- 2 * delta

  if (is.vector(x))
    x <- matrix(x, ncol = length(x))

  d <- ncol(Sigma)

  dec <- tryCatch(Rfast::cholesky(Sigma), error = function(e) e)
  tmp <- backsolve(dec, t(x) - mu, transpose = TRUE)
  rss <- colSums(tmp^2)

  id1 <- which(rss == 0)
  id2 <- which(rss != 0)

  out <- vector("numeric", length(rss))

  out[id1] <- -sum(log(diag(dec))) + log(delta) - (d/2) * log(2 * pi) - log(delta + d)
  out[id2] <- -sum(log(diag(dec))) + log(delta) + (delta / 2 - 1) * log(2) +
    lgamma(0.5 * (delta + d)) + stats::pgamma(rss[id2] / 2, 0.5 * (delta + d), scale = 1, log.p = TRUE) -
    (d/2) * log(pi) - (0.5 * (delta + d)) * log(rss[id2])


  if (!log) out <- exp(out)

  out

}

# Random generation
rmv_slash <- function(n, mu, Sigma, delta) {

  d <- ncol(Sigma)
  A <- Rfast::cholesky(Sigma)
  Z <- matrix(stats::rnorm(n * d), ncol = d)

  mu + Z%*%A / sqrt(stats::rbeta(n, delta, 1))

}

# Distribution function of the squared Mahalanobis distance
maha_slash <- function(q, delta, d) {
  stats::pchisq(q, d) - 2^delta * gamma(delta + d / 2) * stats::pchisq(q, d + 2 * delta) /
    (q^delta * gamma(d / 2))
}

## Multivariate hyperbolic -------------------------------------------------------------------------

# X = sqrt(W) Z,     Z ~ Nd(mu, Sigma) and W ~ GIG(lambda = 1, chi = 1, psi = delta^2)
# Source: https://cran.r-project.org/web/packages/ghyp/vignettes/Generalized_Hyperbolic_Distribution.pdf

# Probability density function
dmv_hyp <- function(x, mu, Sigma, delta, log = FALSE) {

  if (is.vector(x))
    x <- matrix(x, ncol = length(x))

  d <- ncol(Sigma)

  dec <- tryCatch(Rfast::cholesky(Sigma), error = function(e) e)
  tmp <- backsolve(dec, t(x) - mu, transpose = TRUE)
  rss <- colSums(tmp^2)

  out <- -sum(log(diag(dec))) + (d - 1) * log(delta) + log(besselK(delta * sqrt(1 + rss), nu = 1 - d/2)) -
    (d/2) * log(2 * pi) - log(besselK(delta, nu = 1)) -
    (d/2 - 1) * log(delta * sqrt(1 + rss))


  if (!log) out <- exp(out)

  out

}

# Random generation
rmv_hyp <- function(n, mu, Sigma, delta) {

  d <- ncol(Sigma)
  A <- Rfast::cholesky(Sigma)
  Z <- matrix(stats::rnorm(n * d), ncol = d)

  mu + Z%*%A * sqrt(GIGrvg::rgig(n, 1.0, 1.0, delta^2))

}

# Distribution function of the squared Mahalanobis distance
maha_hyp <- function(q, delta, d){


  aux_f <- function(x){

    integrand <- function(u){

      exp((d / 2) * (log(x) - log(2)) - lgamma(d / 2) +
            log(ghyp::pgig(1 / u, 1, 1, delta^2)) + (0.5 * d - 1) * log(u) -0.5 * u * x)

    }

    stats::integrate(integrand, .Machine$double.eps^(1/2), Inf)$val

  }

  apply(matrix(q), 1, aux_f)

}





