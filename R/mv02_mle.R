bcnsm_mle <- function(y, X = NULL,
                      association = c("unstructured", "uniform", "nonassociative"),
                      copula = "gaussian", delta = NULL,
                      margins = "bcno", mu.link = "log",
                      control = control_fit(...), ...)
{

  ### Specific optim parameters
  method <- control$method
  maxit <- control$maxit
  hessian <- control$hessian
  trace <- control$trace

  ### Initial values
  start <- control$start
  mu_inits  <- control$mu_inits
  sigma_inits  <- control$sigma_inits
  lambda_inits <- control$lambda_inits
  nu_inits <- control$nu_inits
  gamma_inits <- control$gamma_inits

  start_id <- !all(is.null(mu_inits), is.null(sigma_inits), is.null(lambda_inits),
                   (nu_inits), is.null(gamma_inits))

  ### Estimation settings
  if (is.vector(y))
    y <- matrix(y, ncol = length(y))

  if (is.data.frame(y))
    y <- as.matrix(y)

  n <- nrow(y)
  d <- ncol(y)

  margins <- as.vector(matrix(margins, 1, d))

  lambda_id <- apply(matrix(margins, ncol = 1), 1, function(x) !grepl("Log-", as.bcs(x)$name))
  nu_id <- apply(matrix(margins, ncol = 1), 1, function(x) as.bcs(x)$extrap)

  control$method <- control$hessian <- control$start <- control$gamma_inits <-
    control$mu_inits <- control$sigma_inits <- control$lambda_inits <- control$nu_inits <- NULL

  ### Covariate matrix
  if (is.null(X)) X <- matrix(1, nrow = n)
  b <- ncol(X)

  ### Copula
  copula <- match.arg(copula, c("gaussian", "t", "slash", "hyp"))

  ### Association matrix
  association <- match.arg(association, c("unstructured", "uniform", "nonassociative"))
  association <- get(association)(d)

  ## Copula specifications
  switch (copula,
          gaussian = {

            dPSI <- function(x, log = FALSE) stats::dnorm(x, log = log)
            qPSI <- function(p) stats::qnorm(p)
            dgf <- function(u, d, log = FALSE){
              out <- -0.5 * u - (d/2) * log(2 * pi)

              if (!log) out <- exp(out)

              out
            }

          },

          t = {

            dPSI <- function(x, log = FALSE) stats::dt(x, delta, log = log)
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
            qPSI <- function(p) distr::q(Whp)(p)
            dgf <- function(u, d, log = FALSE){

              out <- (d - 1) * log(delta) + log(besselK(delta * sqrt(1 + u), nu = 1 - d/2)) -
                (d/2) * log(2 * pi) - log(besselK(delta, nu = 1)) -
                (d/2 - 1) * log(delta * sqrt(1 + u))

              if (!log) out <- exp(out)

              out
            }

          }
  )

  # Parameter indexation -------------------------------------------------------
  par_id <- list(beta = 1 : (d * b),
                 sigma = 1 : d + d * b,
                 lambda = if(any(lambda_id)) 1:(sum(as.numeric(lambda_id))) + d * (b + 1) else NULL,
                 nu = if(any(nu_id)) 1:sum(as.numeric(nu_id)) +
                   (sum(as.numeric(lambda_id))) + d * (b + 1) else NULL,
                 gamma = if(association$npar > 0){
                   1:association$npar + sum(as.numeric(nu_id)) +
                     (sum(as.numeric(lambda_id))) + d * (b + 1)
                 } else {
                   NULL
                 })

  ## Initial values ------------------------------------------------------------
  if (!start_id) {

    if (is.null(gamma_inits)) gamma_inits <- association$start(y)

    beta0 <- matrix(NA, b, d)
    sigma0 <- rep(NA, d)
    lambda0 <- rep(NA, d)
    nu0 <- rep(NA, d)

    for (j in 1:d) {

      if (margins[j] %in% c("bcloi", "lloi", "bcpe", "lpe", "bchp", "lhp", "bcla", "lla")){
        faux <- gamlss.dist::BCPE
      } else if (margins[j] %in% c("bcno", "lno", "bcloii", "lloii")){
        faux <- gamlss.dist::BCCG
      } else{
        faux <- gamlss.dist::BCT
      }

      gamlss_fit <- switch (mu.link,

                            log = {
                              suppressWarnings(gamlss::gamlss(y[,j] ~ X + 0,
                                                              family = faux(mu.link = "log",sigma.link = "identity",
                                                                            nu.link = "identity"),
                                                              trace = FALSE))
                            },

                            identity = {
                              suppressWarnings(gamlss::gamlss(y[,j] ~ X + 0,
                                                              family = faux(mu.link = "identity", sigma.link = "identity",
                                                                            nu.link = "identity"),
                                                              trace = FALSE))
                            }
      )

      beta0[, j] <- stats::coef(gamlss_fit, "mu")
      sigma0[j] <- stats::coef(gamlss_fit, "sigma")
      if (lambda_id[j]) lambda0[j] <- stats::coef(gamlss_fit, "nu")
      if (nu_id[j]) nu0[j] <- min(exp(stats::coef(gamlss_fit, "tau")), 20)

    }

    inits <- c(beta0, sigma0, lambda0[lambda_id], nu0[nu_id], gamma_inits)

  } else {

    inits <- start

  }

  ## Log-likelihood ------------------------------------------------------------
  EPS <- .Machine$double.eps^(1/1.5)

  ll <- function(theta){

    ### Parameter setting
    beta <- matrix(theta[par_id$beta], b, d)

    mu <- stats::make.link(mu.link)$linkinv(X%*%beta)
    sigma <- theta[par_id$sigma]
    lambda <- rep(NA, d)
    lambda[lambda_id] <- theta[par_id$lambda]
    nu <- rep(NA, d)
    nu[nu_id] <- theta[par_id$nu]
    gamma <- theta[par_id$gamma]

    ### Association matrix
    gamma_id <- FALSE
    Gamma <- association$Gamma(gamma)
    dec <- tryCatch(Rfast::cholesky(Gamma), error = function(e) e)
    if (inherits(dec, "error")) dec <- NULL
    if (length(gamma) > 0L) gamma_id <- any(gamma < association$lower | gamma > association$upper)

    ### Out
    if (any(!is.finite(mu)) | any(!is.finite(sigma)) | any(!is.finite(lambda[lambda_id])) |
        any(mu < 0) | any(sigma < 0) | any(nu[nu_id] < 0) | gamma_id | is.null(dec)) {

      -Inf

    }else {

      epsilon <- log_f <- matrix(NA, n, d)
      for(j in 1:d){

        epsilon[, j] <- qPSI(pmin(pmax(get(paste0("p", margins[j]))(q = y[, j],
                                                         mu = mu[, j],
                                                         sigma = sigma[j],
                                                         lambda = lambda[j],
                                                         nu = nu[j]), EPS), 1 - EPS))

        log_f[, j] <- get(paste0("d", margins[j]))(x = y[, j],
                                                   mu = mu[, j],
                                                   sigma = sigma[j],
                                                   lambda = lambda[j],
                                                   nu = nu[j], log = TRUE)
      }

      tmp <- backsolve(dec, t(epsilon), transpose = TRUE)
      rss <- colSums(tmp^2)
      -n * sum(log(diag(dec))) + sum(dgf(rss, d, log = TRUE)) -
        sum(dgf(epsilon^2, 1L, log = TRUE)) + sum(log_f)

    }

  }

  ## Estimates -----------------------------------------------------------------
  opt <- stats::optim(par = inits,
                      fn = ll,
                      method = method,
                      control = control,
                      hessian = FALSE)

  # Assymptotic covariance estimates
  if (hessian) {
    J <- -numDeriv::hessian(ll, opt$par)
    vcov <- try(solve(J), silent = TRUE)
    vcov <- if (unique(grepl("Error", vcov))) matrix(NA, nrow = length(opt$par), ncol = length(opt$par)) else vcov
  }

  if (opt$convergence > 0)
    warning(cat("optimization failed to converge\n"))

  opt$par_id <- par_id
  opt$start <- inits

  beta <- matrix(opt$par[par_id$beta], b, d)
  colnames(beta) <- colnames(y)
  rownames(beta) <- colnames(X)

  sigma <- opt$par[par_id$sigma]
  names(sigma) <- colnames(y)

  lambda <- rep(NA, d)
  names(lambda) <- colnames(y)
  lambda[lambda_id] <- opt$par[par_id$lambda]

  nu <- rep(NA, d)
  names(nu) <- colnames(y)
  nu[nu_id] <- opt$par[par_id$nu]

  # Out list
  out <- list(
    coefficients = beta,
    fitted.values = list(mu = stats::make.link(mu.link)$linkinv(X%*%beta),
                         sigma = sigma, lambda = lambda, nu = nu),
    margins = margins,
    mu.link = mu.link,
    copula = copula,
    delta = delta,
    gamma = if (association$npar > 0L) opt$par[par_id$gamma] else NULL,
    association = association$name,
    logLik = opt$value,
    vcov = if (hessian) vcov else NULL,
    y = y,
    X = X,
    optim_params = opt,
    nobs = n,
    d = d
  )

  out

}
