mr_mle <- function(y, X = NULL, Z = NULL, association = nonassociative(length(y)),
                   mu.link = "log", sigma.link = "log",
                   copula = "gaussian", delta = NULL, margin = "bcno",
                   control = control_fit(...), ...)
{

  ### Specific optim parameters
  method <- control$method
  maxit <- control$maxit
  hessian <- control$hessian

  ### Initial values
  mu_inits  <- control$mu_inits
  sigma_inits  <- control$sigma_inits
  lambda_inits <- control$lambda_inits
  nu_inits <- control$nu_inits
  gamma_inits <- control$gamma_inits

  start_id <- !all(is.null(mu_inits), is.null(sigma_inits), is.null(lambda_inits),
                   (nu_inits), is.null(gamma_inits))

  lambda_id <- !grepl("Log-", as.bcs(margin)$name)
  nu_id <- as.bcs(margin)$extrap

  control$method <- control$hessian <- control$start <- control$gamma_inits <-
    control$mu_inits <- control$sigma_inits <- control$lambda_inits <- control$nu_inits <- NULL

  ### Data setting -------------------------------------------------------------
  n <- length(y)

  if (is.null(X)) X <- matrix(1, nrow = n)
  if (is.null(Z)) Z <- matrix(1, nrow = n)

  k1 <- ncol(X)
  k2 <- ncol(Z)

  ### Copula -------------------------------------------------------------------
  copula <- match.arg(copula, c("gaussian", "t", "slash", "hyp"))

  if (association$name == "non-associative") association$Gamma <- function(gamma = NULL, n) diag(n)

  ### Copula
  copula <- match.arg(copula, c("gaussian", "t", "slash", "hyp"))

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
  d <- 1
  k3 <- 1
  param_id <- list(beta = 1 : (d * k1),
                   kappa = 1 : (d * k2) + d * k1,
                   lambda = if(lambda_id) 1 : (d * k3) + d * (k1 + k2) else NULL,
                   nu = if(any(nu_id)) 1:sum(as.numeric(nu_id)) + d * (k1 + k2 + k3 * as.numeric(lambda_id)) else NULL,
                   gamma = if(association$npar > 0){
                     1:association$npar +
                       sum(as.numeric(nu_id)) + d * (k1 + k2 + k3 * as.numeric(lambda_id))
                   } else {
                     NULL
                   })

  ## Initial values ------------------------------------------------------------
  if (is.null(gamma_inits)) gamma_inits <- association$start(y)

  links <- paste0(mu.link, "_", sigma.link)

  if (margin %in% c("bcpe", "bchp", "bcloi", "lpe", "lhp", "lloi")) {
    faux <- gamlss.dist::BCPE
    faux2 <- function(z, nu) gamlss.dist::pPE(z, 0, 1, nu)
  } else if (margin %in% c("bcno", "lno")) {
    faux <- gamlss.dist::BCCG
    faux2 <- function(z, nu) stats::pnorm(z)
  } else {
    faux <- gamlss.dist::BCT
    faux2 <- function(z, nu) stats::pt(z, nu)
  }


  gamlss_fit <- switch (links,
                        log_log = {
                          suppressWarnings(gamlss::gamlss(y ~ X + 0, ~ Z + 0,
                                                          family = faux(mu.link = "log", sigma.link = "log",
                                                                        nu.link = "identity"),
                                                          trace = FALSE))#, tau.fix = TRUE, tau.start = 4))
                        },

                        log_identity = {
                          suppressWarnings(gamlss::gamlss(y ~ X + 0, ~ Z + 0,
                                                          family = faux(mu.link = "log",sigma.link = "identity",
                                                                        nu.link = "identity"),
                                                          trace = FALSE))#, tau.fix = TRUE, tau.start = 4))
                        },

                        identity_log = {
                          suppressWarnings(gamlss::gamlss(y ~ X + 0, ~ Z + 0,
                                                          family = faux(mu.link = "identity", sigma.link = "log",
                                                                        nu.link = "identity"),
                                                          trace = FALSE))#, tau.fix = TRUE, tau.start = 4))
                        },

                        identity_identity = {
                          suppressWarnings(gamlss::gamlss(y ~ X + 0, ~ Z + 0,
                                                          family = faux(mu.link = "identity", sigma.link = "identity",
                                                                        nu.link = "identity"),
                                                          trace = FALSE))#, tau.fix = TRUE, tau.start = 4))
                        })

  inits <- c(stats::coef(gamlss_fit, "mu"),
             stats::coef(gamlss_fit, "sigma"),
             if(any(lambda_id)) stats::coef(gamlss_fit, "nu") else NULL,
             if(any(nu_id)) min(gamlss_fit$tau.fv[1], 50) else NULL, #max(exp(stats::coef(gamlss_fit, "tau")), 20) else NULL,
             gamma_inits)

  ## Log-likelihood ------------------------------------------------------------
  EPS <- .Machine$double.eps^(1/1.5)

  ll <- function(theta){

    ### Parameter setting ------------------------------------------------------
    beta <- theta[param_id$beta]
    kappa <- theta[param_id$kappa]
    gamma <- if (association$npar > 0) theta[param_id$gamma] else NULL

    # Marginal parameters
    mu <- c(stats::make.link(mu.link)$linkinv(X%*%beta))
    sigma <- c(stats::make.link(sigma.link)$linkinv(Z%*%kappa))
    lambda <- if (lambda_id) theta[param_id$lambda] else NA
    nu <- if (nu_id) theta[param_id$nu] else NA

    ### Association matrix
    Gamma <- association$Gamma(gamma, n)
    dec <- tryCatch(Rfast::cholesky(Gamma), error = function(e) e)
    if (inherits(dec, "error")) dec <- NULL

    ### Out
    if (any(!is.finite(mu)) | any(!is.finite(sigma)) | any(!is.finite(lambda[lambda_id])) |
        any(mu < 0) | any(sigma < 0) | any(nu[nu_id] < 0) | is.null(dec)) {

      -Inf

    }else {

      epsilon <- matrix(qPSI(pmin(pmax(get(paste0("p", margin))(q = y,
                                                         mu = mu,
                                                         sigma = sigma,
                                                         lambda = lambda,
                                                         nu = nu), EPS), 1 - EPS)),
                        nrow = 1)

      log_f <- get(paste0("d", margin))(x = y,
                                        mu = mu,
                                        sigma = sigma,
                                        lambda = lambda,
                                        nu = nu, log = TRUE)


      tmp <- backsolve(dec, t(epsilon), transpose = TRUE)
      rss <- colSums(tmp^2)
      - sum(log(diag(dec))) + dgf(rss, n, log = TRUE) -
        sum(dgf(epsilon^2, 1L, log = TRUE)) + sum(log_f)

    }

  }


  ## Estimates -----------------------------------------------------------------
  opt <- stats::optim(par = inits,
                      fn = ll,
                      method = method,
                      control = control,
                      hessian = hessian)

  opt$par_id <- param_id
  opt$start <- inits

  # Convergence status
  if (opt$convergence > 0)
    warning(cat("optimization failed to converge\n"))


  # Coefficients and marginal parameters
  beta <- opt$par[param_id$beta]
  names(beta) <- colnames(X)
  kappa <- opt$par[param_id$kappa]
  names(kappa) <- colnames(Z)

  mu <- c(stats::make.link(mu.link)$linkinv(X%*%beta))
  sigma <- c(stats::make.link(sigma.link)$linkinv(Z%*%kappa))
  lambda <- if (lambda_id) opt$par[param_id$lambda] else NULL
  nu <- if (nu_id) opt$par[param_id$nu] else NULL

  # Association matrix parameters
  gamma <- if (association$npar) opt$par[param_id$gamma] else NULL

  # Assymptotic covariance estimates
  if (hessian){

    vcov <- try(chol2inv(Rfast::cholesky(-opt$hessian)), silent = TRUE)
    vcov <- if (unique(grepl("Error", vcov))) matrix(NA, nrow = length(opt$par), ncol = length(opt$par)) else vcov

    colnames(vcov) <- rownames(vcov) <- c(colnames(X), colnames(Z),
                                          if (lambda_id) "lambda" else NULL,
                                          if (nu_id) "nu" else NULL,
                                          if (length(gamma)) names(gamma) else NULL)
  }

  out <- list(coefficients = list(beta = beta,
                                 kappa = kappa),
             fitted.values = list(mu = mu,
                                  sigma = sigma,
                                  lambda = lambda,
                                  nu = nu),
             links = list(mu.link = mu.link, sigma.link = sigma.link),
             margin = margin,
             logLik = opt$value,
             vcov = vcov,
             copula = copula,
             delta = delta,
             association = association,
             gamma = gamma,
             nobs = n,
             optim_params = opt,
             y = y, X = X, Z = Z)

  out

}




