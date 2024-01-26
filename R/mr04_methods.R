#' @name marginreg-methods
#' @title Methods for "marginreg" objects
#'
#' @param x,object an object of class \code{marginreg}.
#' @param k numeric, the penalty per parameter to be used; the default
#'     \code{k = 2} is the classical AIC.
#' @param formula a model formula or terms object or an R object.
#' @param ... further arguments passed to or from other methods.
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
NULL


## Model frame
#' @export
#' @rdname marginreg-methods
model.frame.marginreg <- function(formula, ...) {
  formula$terms <- formula$terms$full
  formula$call$formula <- formula$formula <- formula(formula$terms)
  NextMethod()
}

## Model matrix
#' @export
#' @rdname marginreg-methods
model.matrix.marginreg <- function(object, model = c("mu", "sigma"), ...) {
  model <- match.arg(model, c("mu", "sigma"))
  rval <- if(!is.null(object$x[[model]])) object$x[[model]]
  else stats::model.matrix(object$terms[[model]], stats::model.frame(object))
  rval
}

# Parameter estimates
#' @rdname marginreg-methods
#' @export
#' @param model a character indicating which submodel the function should act on. It can be
#'     \code{"mu"} (default) or \code{"sigma"} to specify the regression setings associated
#'     with the \code{mu} or \code{sigma}, respectively.
coef.marginreg <- function(object, model = c("mu", "sigma"), ...) {

  model <- match.arg(model, c("mu", "sigma"))

  switch (model,
          mu = object$coefficients$beta,
          sigma = object$coefficients$kappa)

}

#  Variance-covariance matrix
#' @rdname marginreg-methods
#' @param parameter character; specifies which submatrix of the asymptotic covariance matrix of the
#'     maximum likelihood estimators should be returned. The options are \code{"full"} (default),
#'     \code{"mu"}, \code{"sigma"}, \code{"lambda"}, \code{"nu"}, and \code{"gamma"}.
#' @export
vcov.marginreg <- function(object, parameter = c("full", "mu", "sigma", "lambda", "nu", "gamma"), ...) {

  parameter <- match.arg(parameter, c("full", "mu", "sigma", "lambda", "nu", "gamma"))
  covm <- object$vcov

  margin <- object$margin
  lambda_id <- !grepl("Log-", as.bcs(margin)$name)
  nu_id <- as.bcs(margin)$extrap

  par_id <- object$optim_params$par_id

  switch(parameter,
         "full" = covm,
         "mu" = covm[par_id$beta, par_id$beta],
         "sigma" = covm[par_id$kappa, par_id$kappa],
         "lambda" = covm[par_id$lambda, par_id$lambda],
         "nu" = covm[par_id$nu, par_id$nu],
         "gamma" = covm[par_id$gamma, par_id$gamma]
  )

}


# Log-likelihood
#' @rdname marginreg-methods
#' @export
logLik.marginreg <- function(object, ...) {
  structure(object$logLik,
            df = length(object$optim$par),
            class = "logLik")
}


# AIC
#' @export
#' @rdname marginreg-methods
AIC.marginreg <- function(object, ..., k = 2) {

  npar <- length(object$optim$par)
  AIC <- - 2 * object$logLik + k * (npar + !is.null(object$delta))

  class(AIC) <- "AIC"
  return(AIC)
}

# Residuals
#' @name residuals.marginreg
#' @title Extract Model Residuals
#'
#' @param object a \code{"marginreg"} object.
#' @param type character; it specifies which residual should be extracted.
#'     The available arguments are \code{"conditional"} (default) and \code{"quantile"}.
#' @param ... further arguments passed to or from other methods.
#'
#' @export
#'
residuals.marginreg <- function(object, type = c("conditional", "quantile"), ...)
{

  type <- match.arg(type, c("conditional", "quantile"))

  copula <- object$copula
  delta <- object$delta
  association <- object$association
  gamma <- object$gamma

  ## Univariate distribution and quantile function of the corresponding NSM distribution
  if (copula == "slash") {
    Wsl <- distr::AbscontDistribution(
      d = function(x) dslash(x, nu = delta),
      Symmetry = distr::SphericalSymmetry(0)
    )
    pPSI <- function(q) distr::p(Wsl)(q)
    qPSI <- function(p) distr::q(Wsl)(p)
  } else if (copula == "gaussian") {
    pPSI <- function(q) stats::pnorm(q)
    qPSI <- function(p) stats::qnorm(p)
  } else if (copula == "t") {
    pPSI <- function(q) stats::pt(q, delta)
    qPSI <- function(p) stats::qt(p, delta)
  } else if (copula == "hyp"){
    Whp <- distr::AbscontDistribution(
      d = function(x) dhyp(x, nu = delta),
      Symmetry = distr::SphericalSymmetry(0)
    )
    pPSI <- function(q) distr::p(Whp)(q)
    qPSI <- function(p) distr::q(Whp)(p)
  }

  n <- object$nobs

  margin <- object$margin

  lambda_id <- !grepl("Log-", as.bcs(margin)$name)
  nu_id <- as.bcs(margin)$extrap

  pBCS <- get(paste0("p", margin))

  mu <- object$fitted.values$mu
  sigma <- object$fitted.values$sigma
  lambda <- if (lambda_id) object$fitted.values$lambda else NA
  nu <- if (nu_id) object$fitted.values$nu else NA

  y <- if(is.null(object$y)) stats::model.response(stats::model.frame(object)) else object$y

  EPS <- .Machine$double.eps^(1/1.5)
  epsilon <- qPSI(pmin(pmax(pBCS(q = y, mu = mu, sigma = sigma, lambda = lambda, nu = nu), EPS), 1 - EPS))

  # A^{-T}%*%epsilon
  Gamma <- association$Gamma(gamma, n)
  A <- Rfast::cholesky(Gamma)
  epsilon_star <- c(as.numeric(solve(t(A))%*%epsilon))

  if (type == "conditional") {

    if(copula == "gaussian"){

      residuals <- epsilon_star

    }else if(copula == "t"){

      Mj <- function(j)
      {
        g  <- sum(epsilon_star[1:(j - 1)]^2)
        s1 <- sqrt((delta + j - 1) / (delta + g))
        stats::pt(s1 * epsilon_star[j], df = delta + j - 1)
      }

      residuals <- c(stats::qnorm(c(stats::pt(epsilon_star[1], delta), apply(matrix(2:n), 1, Mj))))

    }else if(copula == "slash"){

      Mj <- function(j){

        a <- t(epsilon_star[1:(j-1)])%*%epsilon_star[1:(j-1)]

        integrand <- function(x){

          aux_f <- function(u){

            exp( (0.5 * (delta + j - 1)) * log(a) - log(pi)/2 -
                   lgamma(0.5 * (delta + j - 1)) -
                   stats::pgamma(0.5 * a, 0.5 * (delta + j - 1), scale = 1, log.p = TRUE) +
                   lgamma(0.5 * (delta + j)) + stats::pgamma(0.5 * (a + u^2), 0.5 * (delta + j), scale = 1, log.p = TRUE) -
                   (0.5 * (delta + j)) * log(a + u^2) )

          }

          apply(matrix(x), 1, aux_f)

        }

        stats::integrate(integrand, -Inf, epsilon_star[j])$val
      }

      residuals <- c(stats::qnorm(c(pPSI(epsilon_star[1]), apply(matrix(2:n), 1, Mj))))

    }else if(copula == "hyp"){

      Mj <- function(j){

        a <- t(epsilon_star[1:(j-1)])%*%epsilon_star[1:(j-1)]

        integrand <- function(x){

          aux_f <- function(u){

            exp( log(delta) - 0.5 * log(2 * pi) +
                   log(besselK(delta * sqrt(a + u^2), 1 - j/2)) -
                   log(besselK(delta * sqrt(a), 1 - (j-1)/2)) +
                   (0.5 * (j - 1)  - 1) * log(delta * sqrt(a)) -
                   (0.5 * j  - 1) * log(delta * sqrt(a + u^2)))

          }

          apply(matrix(x), 1, aux_f)

        }

        stats::integrate(integrand, -Inf, epsilon_star[j])$val
      }

      residuals <- c(stats::qnorm(c(pPSI(epsilon_star[1]), apply(matrix(2:n), 1, Mj))))

    }

  } else if (type == "quantile") {

    residuals <- c(stats::qnorm(pBCS(q = y, mu = c(mu), sigma = c(sigma), lambda = c(lambda), nu = nu)))

  }

  residuals

}


# Print
#' @rdname marginreg-methods
#' @export
print.marginreg <- function(x, ...)
{

  n <- x$nobs

  ## Link function
  mu.link <- x$link$mu.link
  sigma.link <- x$link$sigma.link

  margin <- as.bcs(x$margin)
  lambda_id <- !grepl("Log-", margin$name)
  nu_id <- margin$extrap
  gamma_id <- x$association$name != "non-associative"

  copula_name <- if (x$copula == "gaussian") "Gaussian" else x$copula

  if (requireNamespace("crayon", quietly = TRUE)) {

    cat(crayon::cyan(margin$name, "marginal regression with",
                     if (is.null(x$delta)) copula_name else paste0(tolower(copula_name), "(", round(x$delta, 2), ")"),
                     "copula\n"))

    cat(crayon::cyan("\nCall:\n"))
    print(x$call)
    if (x$optim$convergence) {
      cat("\nmodel did not converge\n")
    } else {
      cat(crayon::cyan("\nMu coefficients with", mu.link, "link function:\n"))
      print(round(x$coefficients$beta, 4))

      cat(crayon::cyan("\nSigma coefficients with", sigma.link, "link function:\n"))
      print(round(x$coefficients$kappa, 4))

      if (lambda_id | nu_id) {

        tab <- round(cbind(lambda = c(x$fitted.values$lambda),
                     nu = c(x$fitted.values$nu)), 4)
        colnames(tab) <- c(if(any(lambda_id)) "lambda", if(any(nu_id)) "nu")
        rownames(tab) <- c("")
        cat(crayon::cyan("\nFurther marginal parameters:\n"))
        print(tab, na.print = "-")

      }

      if (gamma_id) {

        if (grepl("ARMA", x$association$name)) {
          cat(crayon::cyan("\nARMA coefficients:\n"))
        } else if(grepl("spatial", x$association$name)) {
          cat(crayon::cyan("\nSpatial coefficient:\n"))
        } else {
          cat(crayon::cyan("\nAssociation coefficients:\n"))
        }

        ngamma <- length(x$gamma)
        npar <- length(x$optim$par)
        tab_gamma <- round(matrix(c(x$gamma), byrow = TRUE, nrow = 1), 4)
        colnames(tab_gamma) <- names(x$gamma)
        rownames(tab_gamma) <- c("")
        print(tab_gamma)

      }

    }


    cat(
      crayon::cyan("\nAssociation:"), tolower(x$association$name),
      crayon::cyan("\nlogLik:"), x$logLik, "|",
      crayon::cyan("AIC:"), stats::AIC(x), "|",
      crayon::cyan("BIC:"), stats::AIC(x, k = log(n)), "\n"
    )

  } else {

    cat(margin$name, "marginal regression with",
        if (is.null(x$delta)) copula_name else paste0(tolower(copula_name), "(", round(x$delta, 2), ")"),
        "copula\n")

    cat("\nCall:\n")
    print(x$call)
    if (x$optim$convergence) {
      cat("\nmodel did not converge\n")
    } else {
      cat("\nMu coefficients with", mu.link, "link function:\n")
      print(round(x$coefficients$beta, 4))

      cat("\nSigma coefficients with", sigma.link, "link function:\n")
      print(round(x$coefficients$kappa, 4))

      if (lambda_id | nu_id) {

        tab <- cbind(lambda = c(x$fitted.values$lambda, sqrt(stats::vcov(x, "lambda"))),
                     nu = c(x$fitted.values$nu, sqrt(stats::vcov(x, "nu"))))
        colnames(tab) <- c(if(any(lambda_id)) "lambda", if(any(nu_id)) "nu")
        rownames(tab) <- c("est.", "s.e.")
        cat("\nFurther marginal parameters:\n")
        print(tab, na.print = "-")

      }

      if (gamma_id) {

        if (grepl("ARMA", x$association$name)) {
          cat("\nARMA coefficients:\n")
        } else if(grepl("spatial", x$association$name)) {
          cat("\nSpatial coefficient:\n")
        } else {
          cat("\nAssociation coefficients:\n")
        }

        ngamma <- length(x$gamma)
        npar <- length(x$optim$par)
        gamma_vcov <- stats::vcov(x, "gamma")

        se_gamma <- if (ngamma == 1) sqrt(gamma_vcov) else sqrt(diag(gamma_vcov))
        tab_gamma <- round(matrix(c(x$gamma, se_gamma), byrow = TRUE, nrow = 2), 4)
        colnames(tab_gamma) <- names(x$gamma)
        rownames(tab_gamma) <- c("est", "s.e.")
        print(tab_gamma)

      }

    }

    cat(
      "\nAssociation:", tolower(x$association$name),
      "\nlogLik:", x$logLik, "|",
      "AIC:", stats::AIC(x), "|",
      "BIC:", stats::AIC(x, k = log(n)), "\n"
    )

  }

  invisible(x)
}


# Summary
#' @rdname marginreg-methods
#' @export
summary.marginreg <- function(object, ...)
{

  n <- object$nobs

  ## Link function
  mu.link <- object$link$mu.link
  sigma.link <- object$link$sigma.link

  margin <- as.bcs(object$margin)
  lambda_id <- !grepl("Log-", margin$name)
  nu_id <- margin$extrap
  gamma_id <- object$association$name != "non-associative"

  # Summary for the conditional quantile residuals
  res <- stats::residuals(object)
  skewness <- mean((res - mean(res))^3) / (stats::sd(res)^3)
  kurtosis <- mean((res - mean(res))^4) / (stats::sd(res)^4)
  tab_residuals <- round(cbind(mean(res), stats::sd(res), skewness, kurtosis), 4)
  colnames(tab_residuals) <- c("Mean", "Std. dev.", "Skewness", "Kurtosis")
  rownames(tab_residuals) <- " "

  # Summary for mu
  est.mu <- stats::coef(object, "mu")
  se.mu <- sqrt(diag(stats::vcov(object, "mu")))
  zval.mu <- est.mu/se.mu
  pval.mu <- 2 * stats::pnorm(abs(zval.mu), lower.tail = FALSE)

  tab_mu <- round(cbind(Estimate = est.mu, `Std. Error` = se.mu,
                  `z value` = zval.mu, `Pr(>|z|)` = pval.mu), 4)

  # Summary for sigma
  est.sigma <- stats::coef(object, "sigma")
  se.sigma <- sqrt(diag(stats::vcov(object, "sigma")))
  zval.sigma <- est.sigma/se.sigma
  pval.sigma <- 2 * stats::pnorm(abs(zval.sigma), lower.tail = FALSE)

  tab_sigma <- round(cbind(Estimate = est.sigma, `Std. Error` = se.sigma,
                     `z value` = zval.sigma, `Pr(>|z|)` = pval.sigma), 4)

  # Summary for lambda or nu
  if (lambda_id | nu_id) {

    tab_lambda_nu <- round(cbind(if (lambda_id) lambda = c(object$fitted.values$lambda,
                                                     sqrt(stats::vcov(object, "lambda"))),
                           if (nu_id) nu = c(object$fitted.values$nu,
                                             sqrt(stats::vcov(object, "nu")))), 4)
    colnames(tab_lambda_nu) <- c(if (lambda_id) "lambda", if(nu_id) "nu")
    rownames(tab_lambda_nu) <- c("est.", "s.e.")
    cat("\nFurther marginal parameters:\n")
    print(tab_lambda_nu, na.print = "-")

  } else {

    tab_lambda_nu <- NULL

  }

  # Summary for gamma
  if (gamma_id) {

    ngamma <- length(object$gamma)
    npar <- length(object$optim$par)
    gamma_vcov <- object$vcov[(npar - ngamma + 1):npar, (npar - ngamma + 1):npar]

    se_gamma <- if (ngamma == 1) sqrt(gamma_vcov) else sqrt(diag(gamma_vcov))
    tab_gamma <- round(matrix(c(object$gamma, se_gamma), byrow = TRUE, nrow = 2), 4)
    colnames(tab_gamma) <- names(object$gamma)
    rownames(tab_gamma) <- c("est.", "s.e.")

  } else {

    tab_gamma <- NULL

  }

  out <- object
  out$tab_residuals <- tab_residuals
  out$tab_mu <- tab_mu
  out$tab_sigma <- tab_sigma
  out$tab_lambda_nu <- tab_lambda_nu
  out$tab_gamma <- tab_gamma
  out$AIC <- stats::AIC(object)
  out$BIC <- stats::AIC(object, k = log(n))

  class(out) <- "summary.marginreg"
  out

}

# Print summary
#' @rdname marginreg-methods
#' @export
print.summary.marginreg <- function(x, ...)
{

  n <- x$nobs

  ## Link function
  mu.link <- x$link$mu.link
  sigma.link <- x$link$sigma.link

  margin <- as.bcs(x$margin)
  lambda_id <- !grepl("Log-", margin$name)
  nu_id <- margin$extrap
  gamma_id <- x$association$name != "non-associative"

  copula_name <- if (x$copula == "gaussian") "Gaussian" else x$copula

  cat(crayon::cyan(margin$name, "marginal regression with",
                   if (is.null(x$delta)) copula_name else paste0(tolower(copula_name), "(", round(x$delta, 2), ")"),
                   "copula\n"))

  cat(crayon::cyan("\nCall:\n\n"))
  print(x$call)

  cat(crayon::cyan("\nConditional quantile residuals:\n\n"))
  print(x$tab_residuals)

  if (x$optim$convergence) {

    cat("\nmodel did not converge\n")

  } else {

    cat(crayon::cyan("\nMu coefficients with", mu.link, "link function:\n\n"))
    stats::printCoefmat(x$tab_mu)

    cat(crayon::cyan("\nSigma coefficients with", sigma.link, "link function:\n\n"))
    stats::printCoefmat(x$tab_sigma)

    if (lambda_id | nu_id) {
      cat(crayon::cyan("\nFurther marginal parameters:\n\n"))
      stats::printCoefmat(x$tab_lambda_nu)
    }

    if (gamma_id) {

      if (grepl("ARMA", x$association$name)) {
        cat(crayon::cyan("\nARMA coefficients:\n\n"))
      } else if(grepl("spatial", x$association$name)) {
        cat(crayon::cyan("\nSpatial coefficient:\n"))
      } else {
        cat(crayon::cyan("\nAssociation coefficients:\n"))
      }

      print(x$tab_gamma)
    }


  }

  cat(
    crayon::cyan("\nAssociation:"), tolower(x$association$name),
    crayon::cyan("\nlogLik:"), x$logLik, "|",
    crayon::cyan("AIC:"), x$AIC, "|",
    crayon::cyan("BIC:"), x$BIC, "\n"
  )

  invisible(x)
}


# Plot
#' @name plot.marginreg
#' @title Diagnostic plots for the BCNSM marginal regression models
#'
#' @param x an object of class \code{"marginreg"}.
#' @param which numeric. If a subset of the plots is required, specify a subset
#'     of the numbers \\code{1:6}. See details below.
#' @param type character; it specifies which residual should be extracted.
#'     The available arguments are \code{"conditional"} (default), \code{"quantile"}, and
#'     \code{"epsilon"} (standardized epsilon's).
#' @param ask logical. If \code{TRUE}, the user is asked before each plot.
#' @param ... further arguments passed to or from other methods.
#'
#' @details The \code{which} argument can be used to select a subset of currently six supported
#'     types of displays. The displays are:
#'
#' \describe{
#'   \item{\code{which = 1}:}{residuals vs fitted \code{mu} values.}
#'   \item{\code{which = 2}:}{residuals vs indices of observations.}
#'   \item{\code{which = 3}:}{residuals density vs standard normal density.}
#'   \item{\code{which = 4}:}{normal probability plot of the residuals.}
#'   \item{\code{which = 5}:}{autocorrelation plot of the residuals.}
#'   \item{\code{which = 6}:}{partial autocorrelation plot of the residuals.}
#'  }
#'
#' @export
#'
plot.marginreg <- function(x, which = 1:4, type = c("conditional", "quantile"),
                      ask = prod(graphics::par("mfcol")) < length(which) &&
                        grDevices::dev.interactive(), ...)
{

  if(!is.numeric(which) || any(which < 1) || any(which > 6))
    stop("`which' must be in 1:6")

  type <- match.arg(type, c("conditional", "quantile"))

  ## Reading
  res <- stats::residuals(x, type)
  n <- x$nobs

  ## Legends
  types <- c("conditional", "quantile")
  Types <- c("Conditional quantile residuals", "Quantile residuals")
  Type <- Types[type == types]

  ## Graphical parameters setting
  if (ask) {
    op <- graphics::par(ask = TRUE)
    on.exit(graphics::par(op))
  }

  ## Plots to shown
  show <- rep(FALSE, 6)
  show[which] <- TRUE

  ## Residuals versus Fitted values
  if (show[1]){
    mu <- stats::fitted.values(x)$mu
    print(ggplot2::ggplot(data.frame(x = mu, y = res), ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_abline(intercept = c(0, stats::qnorm(c(0.025, 1 - 0.025))),
                  slope = 0, lty = c(1, 2, 2), lwd = c(1, 0.5, 0.5), col = "#56B1F7") +
      ggplot2::geom_point(pch = "+", cex = 3) +
      ggplot2::labs(x = "Fitted values", y = Type))

  }

  ## Residuals versus index observation
  if (show[2]){

    print(ggplot2::ggplot(data.frame(x = 1:n, y = res), ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_abline(intercept = c(0, stats::qnorm(c(0.025, 1 - 0.025))),
                           slope = 0, lty = c(1, 2, 2), lwd = c(1, 0.5, 0.5), col = "#56B1F7") +
      ggplot2::geom_point(pch = "+", cex = 3) +
      ggplot2::labs(x = "Fitted values", y = Type))

  }

  ## Normal density comparison
  if(show[3]) {

    print(ggplot2::ggplot(data = data.frame(x = res), ggplot2::aes(x = x)) +
      ggplot2::geom_density() +
      ggplot2::geom_rug() +
      ggplot2::geom_function(fun = stats::dnorm, col = "#56B1F7") +
      ggplot2::xlim(-4, 4) +
      ggplot2::labs(x = Type, y = "Density"))

  }

  ## Normal probability plot
  if(show[4]) {

    Pi <- (1:n - 0.5) / n
    zi <- stats::qnorm(Pi)

    mu_hat <- stats::median(res)
    sigma_hat <- stats::IQR(res) / 1.349

    se <- sigma_hat * sqrt(Pi * (1 - Pi) / n) / stats::dnorm(zi)
    rq_hat <- mu_hat + sigma_hat * zi

    positions <- data.frame(x = c(zi, rev(zi)),
                            y = c(rq_hat - stats::qnorm(1 - 0.05 / 2) * se,
                                  rev(rq_hat + stats::qnorm(1 - 0.05 / 2) * se)))

    print(ggplot2::ggplot(data.frame(zi = zi, res = sort(res))) +
      ggplot2::geom_polygon(data = positions, ggplot2::aes(x = x, y = y), fill = "#cceff1") +
      ggplot2::geom_point(ggplot2::aes(x = zi, y = res), pch = "+", cex = 4) +
      ggplot2::labs(x = "Normal quantiles", y = "Residual quantiles") +
      ggplot2::geom_qq_line(ggplot2::aes(sample = res), col = "#56B1F7", lwd = 1))

  }

  ## ACF of residuals
  if(show[5]) {

    print(forecast::ggAcf(res, lag.max = min(n - 20, 50)) +
      ggplot2::labs(title = "") + ggplot2::ylim(-1, 1))

  }

  ## PACF of residuals
  if(show[6]) {

    print(forecast::ggPacf(res, lag.max = min(n - 20, 50)) +
      ggplot2::labs(title = "") + ggplot2::ylim(-1, 1))

  }

}

globalVariables("y")

#' Predict Method for BCNSM Marginal Regression Fits
#'
#' Obtains predictions from a fitted BCNSM marginal regression object.
#'
#' @param object an \code{"marginreg"} object.
#' @param newdata optionally, a data frame in which to look for variables
#'     with which to predict. If omitted, the fitted quantiles are returned.
#' @param at the order of the quantile to be predicted. The default is to predict the median,
#'     that is, \code{at = 0.5}.
#' @param na.action function determining what should be done with missing
#'     values in \code{newdata}. The default is to predict \code{NA}.
#' @param h number of steps forward for time series forecasting.
#' @param dist a distance matrix used specifically for predictions of spatial data at locations not
#'     originally observed in the data set. Each row of the matrix is assumed to contain the distances
#'     from the original observations to an unobserved location.
#' @param ...  arguments passed to or from other methods.
#'
#' @return A vector of predictions.
#' @export
#'
predict.marginreg <- function(object, newdata = NULL, at = 0.5, na.action = stats::na.pass, h, dist, ...)
{

  if (at <= 0 | at >= 1)
    stop("at must be a scalar on the unit interval (0, 1)")

  n <- object$nobs

  mu <- object$fitted.values$mu
  sigma <- object$fitted.values$sigma
  lambda <- object$fitted.values$lambda
  nu <- object$fitted.values$nu

  margin <- object$margin

  qbcs <- get(paste0("q", margin))
  pbcs <- get(paste0("p", margin))

  if(missing(newdata)) {

    return(qbcs(p = rep(at, n), mu = mu, sigma = sigma, lambda = lambda, nu = nu))

  } else {

    if (copula != "gaussian")
      stop("\npredictions for newdata is currently available only for the Gaussian copula\n")

    copula <- object$copula
    delta <- object$delta
    association <- object$association
    gamma <- object$gamma

    ## Univariate distribution and quantile function of the corresponding NSM distribution
    if (copula == "slash") {
      Wsl <- distr::AbscontDistribution(
        d = function(x) dslash(x, nu = delta),
        Symmetry = distr::SphericalSymmetry(0)
      )
      qPSI <- function(p) distr::q(Wsl)(p)
    } else if (copula == "gaussian") {
      qPSI <- function(p) stats::qnorm(p)
    } else if (copula == "t") {
      qPSI <- function(p) stats::qt(p, delta)
    } else if (copula == "hyp"){
      Whp <- distr::AbscontDistribution(
        d = function(x) dhyp(x, nu = delta),
        Symmetry = distr::SphericalSymmetry(0)
      )
      qPSI <- function(p) distr::q(Whp)(p)
    }

    EPS <- .Machine$double.eps^(1/1.5)
    y <- if(is.null(object$y)) stats::model.response(stats::model.frame(object)) else object$y
    epsilon <- qPSI(pmin(pmax(pbcs(q = y, mu = mu, sigma = sigma, lambda = lambda, nu = nu), EPS), 1 - EPS))

    mf <- stats::model.frame(stats::delete.response(object$terms[["full"]]),
                             newdata, na.action = na.action)
    newdata <- newdata[rownames(mf), , drop = FALSE]

    X_new <- stats::model.matrix(stats::delete.response(object$terms$mu), mf)
    Z_new <- stats::model.matrix(stats::delete.response(object$terms$sigma), mf)

    mu.link <- object$links$mu.link
    sigma.link <- object$links$sigma.link

    beta <- object$coefficients$beta
    kappa <- object$coefficients$kappa

    mu_new <- stats::make.link(mu.link)$linkinv(X_new%*%beta)
    sigma_new <- stats::make.link(sigma.link)$linkinv(Z_new%*%kappa)
    len <- nrow(X_new)

    if (grepl("ARMA", object$association$name)) {

      p <- sum(grepl("ar", names(gamma)))
      q <- sum(grepl("ma", names(gamma)))
      m <- max(p, q + 1)

      PHI <- if(p >= 1L) gamma[1:p] else NULL
      THETA <- if(q >= 1L) gamma[1:q + p] else NULL

      process <- paste0("arma", p, q)
      sigma2 <- switch(process,
                       "arma01" = 1 / (1 + THETA[1]^2),
                       "arma02" = 1 / (1 + THETA[1]^2 + THETA[2]^2),
                       "arma10" = 1 - PHI[1]^2,
                       "arma11" = (1 - PHI[1]^2) / (1 + 2 * PHI[1] * THETA[1] + THETA[1]^2),
                       "arma12" = (1 - PHI[1]^2) / (1 + THETA[1]^2 + THETA[2]^2 +
                                                      PHI[1] * THETA[1] +
                                                      PHI[1] * THETA[1] * THETA[2] +
                                                      PHI[1]^2 * THETA[2]),
                       "arma20" = (1 + PHI[2]) * (1 - PHI[2] - PHI[1]) * (1 - PHI[2] + PHI[1]) / (1 - PHI[2]),
                       "arma21" = (1 + PHI[2]) * (1 - PHI[2] - PHI[1]) * (1 - PHI[2] + PHI[1]) /
                         (1 + THETA[1]^2 + 2 * PHI[1] * THETA[1] - PHI[2] * THETA[1]^2 - PHI[2]),
                       "arma22" = (1 + PHI[2]) * (1 - PHI[2] - PHI[1]) * (1 - PHI[2] + PHI[1]) /
                         (2 * PHI[1]^2 * THETA[2] + 2 * PHI[1] * THETA[1] * (1 + THETA[2]) +
                            (1 - PHI[2]) * (1 + THETA[1]^2 + THETA[2]^2 + 2 * PHI[2] * THETA[2])))

      # State space model representation
      SSM <- dlm::dlmModARMA(ar = PHI, ma = THETA, sigma2 = sigma2, dV = 0L)
      filter <- dlm::dlmFilter(epsilon, SSM)
      forecast <- dlm::dlmForecast(filter, nAhead = h)
      mh <- c(forecast$f)
      sh <- sqrt(unlist(forecast$Q))

      pred <- vector("numeric", h)
      for(i in 1:h){

        w <- pmin(pmax(stats::pnorm(mh[i] + stats::qnorm(at) * sh[i]),
                       .Machine$double.eps^(1/2)), 1 - .Machine$double.eps^(1/2))
        pred[i] <- qbcs(p = w,
                        mu = mu_new[i],
                        sigma = sigma_new[i],
                        lambda = lambda,
                        nu = nu)

      }


    } else if (grepl("spatial", object$association$name)) {

      pred <- vector("numeric", length = len)
      for (i in 1:len){

        Gamma11 <- association$Gamma(gamma, n)
        Gamma12 <- matrix(.matern(dist[i, ]))
        Gamma21 <- t(Gamma12)

        Gamma11_inv <- chol2inv(Rfast::cholesky(as.matrix(Gamma11)))

        mu_star <- Gamma21%*%Gamma11_inv%*%epsilon
        sigma_star <- 1 - Gamma21%*%Gamma11_inv%*%Gamma12

        w <- pmin(pmax(stats::pnorm(mu_star + stats::qnorm(at) * sigma_star),
                       .Machine$double.eps^(1/2)), 1 - .Machine$double.eps^(1/2))

        pred[i] <- qbcs(p = w, mu = mu_new[i], sigma = sigma_new[i], lambda = lambda, nu = nu)

      }

    } else {

      pred <- qbcs(p = rep(at, length(mu_new)), mu = mu_new, sigma = sigma_new, lambda = lambda, nu = nu)

    }


    return(pred)

  }
}
