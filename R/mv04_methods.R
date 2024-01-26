#' @name multireg-methods
#' @title Methods for "multireg" objects
#' @param x,object an object of class \code{"multireg"}.
#' @param formula a simbolic description of the model to be fitted, of type
#'     \code{cbind(y1, y2, ..., yd) ~ x1 + x2 + ... + xb}.
#' @param k numeric, the penalty per parameter to be used; the default
#'     \code{k = 2} is the classical AIC.
#' @param ... further arguments passed to or from other methods.
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
NULL

## Model frame
#' @export
#' @rdname multireg-methods
model.frame.multireg <- function(formula, ...) {
  formula$terms <- formula$terms
  formula$call$formula <- formula$formula <- formula(formula$terms)
  NextMethod()
}

## Model matrix
#' @export
#' @rdname multireg-methods
model.matrix.multireg <- function(object, ...) {
  rval <- if(!is.null(object$x)) object$x
  else stats::model.matrix(object$terms, stats::model.frame(object))
  return(rval)
}


#  Variance-covariance matrix
#' @rdname multireg-methods
#' @param parameter character; specifies which submatrix of the asymptotic covariance matrix of the
#'     maximum likelihood estimators should be returned. The options are \code{"full"} (default),
#'     \code{"beta"}, \code{"sigma"}, \code{"lambda"}, \code{"nu"}, and \code{"gamma"}.
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @export
vcov.multireg <- function(object, parameter = c("full", "beta", "sigma", "lambda", "nu", "gamma"), ...) {

  parameter <- match.arg(parameter, c("full", "beta", "sigma", "lambda", "nu", "gamma"))
  covm <- object$vcov

  margins <- object$margins
  lambda_id <- apply(matrix(margins, ncol = 1), 1, function(x) !grepl("Log-", as.bcs(x)$name))
  nu_id <- apply(matrix(margins, ncol = 1), 1, function(x) as.bcs(x)$extrap)

  y <- if (is.null(object$y)) stats::model.response(stats::model.frame(object))  else object$y
  X <- stats::model.matrix(object)

  d <- object$d
  b <- ncol(X)

  names_beta <- paste0(rep(colnames(y), each = ncol(X)), "_", rep(colnames(X), d))

  if (object$association == "uniform") {
    names_gamma <- "gamma"
  } else if (object$association == "unstructured") {
    names_gamma <- vector()
    id <- matrix(which(upper.tri(diag(d)), arr.ind = TRUE)[order(which(upper.tri(diag(d)), arr.ind = TRUE)[, 1]), ], ncol = 2)
    for (i in 1:nrow(id)) {
      names_gamma[i] <- paste0("gamma", id[i, 1], id[i, 2])
    }
  }


  colnames(covm) <- rownames(covm) <- c(
    names_beta,
    paste0("sigma", 1:d),
    if (any(lambda_id)) paste0("lambda", (1:d)[lambda_id]) else NULL,
    if (any(nu_id)) paste0("nu", (1:d)[nu_id]) else NULL,
    names_gamma
  )

  par_id <- object$optim_params$par_id

  covm_beta <- list()
  for (j in 1:d) {
    covm_beta[[colnames(y)[j]]] <- covm[1:b + 2 * (j - 1), 1:b + 2 * (j - 1)]
  }

  switch(parameter,
         "full" = covm,
         "beta" = covm_beta,
         "sigma" = covm[par_id$sigma, par_id$sigma],
         "lambda" = covm[par_id$lambda, par_id$lambda],
         "nu" = covm[par_id$nu, par_id$nu],
         "gamma" = covm[par_id$gamma, par_id$gamma]
  )
}

# Log-likelihood
#' @rdname multireg-methods
#' @export
logLik.multireg <- function(object, ...) {
  structure(object$logLik,
            df = length(object$optim_params$par) + as.numeric(!is.null(object$delta)),
            class = "logLik")
}

# AIC
#' @export
#' @rdname multireg-methods
AIC.multireg <- function(object, ..., k = 2) {
  npars <- length(object$optim_params$par) + as.numeric(!is.null(object$delta))
  AIC <- -2 * object$logLik + k * npars

  class(AIC) <- "AIC"
  return(AIC)
}

# Residuals
#' @name residuals.multireg
#' @title Extract Model Residuals for a Fitted BCNSM Regression
#'
#' @param object an \code{"multireg"} object.
#' @param type character; specifies which residual should be extracted. The available arguments are
#'     \code{"quantile"} (default), \code{"mahalanobis"}, and \code{"epsilon"}.
#' @param ... further arguments passed to or from other methods.
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @export
#'
residuals.multireg <- function(object, type = c("quantile", "marginal", "epsilon"), ...) {

  ## raw response residuals and desired type
  type <- match.arg(type, c("quantile", "marginal", "epsilon"))

  y <- if(is.null(object$y)){
    stats::model.response(stats::model.frame(object))
  } else{
    object$y
  }

  n <- object$nobs
  d <- object$d
  margins <- object$margins
  copula <- object$copula
  delta <- object$delta

  # Univariate marginal quantile function of the NSM distribution
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

  } else {
    Whp <- distr::AbscontDistribution(
      d = function(x) dhyp(x, nu = delta),
      Symmetry = distr::SphericalSymmetry(0)
    )
    qPSI <- function(p) distr::q(Whp)(p)
  }

  mu <- object$fitted.values$mu
  sigma <- object$fitted.values$sigma
  lambda <- object$fitted.values$lambda
  nu <- object$fitted.values$nu

  EPS <- .Machine$double.eps

  epsilon <- rm <- matrix(NA, n, d)
  for(j in 1:d){

    epsilon[, j] <- qPSI(pmin(pmax(get(paste0("p", margins[j]))(as.matrix(y[, j]), mu = mu[, j],
                                                                sigma = sigma[j], lambda = lambda[j],
                                                                nu = nu[j]), EPS), 1 - EPS))

    rm[, j] <- stats::qnorm(pmin(pmax(get(paste0("p", margins[j]))(as.matrix(y[, j]), mu = mu[, j],
                                                            sigma = sigma[j], lambda = lambda[j],
                                                            nu = nu[j]), EPS), 1 - EPS))

  }

  colnames(epsilon) <- colnames(rm) <- colnames(y)

  Gamma <- get(object$association)(d)$Gamma(object$gamma)

  # Squared Mahalanobis distance
  mahalanobis <- mahalanobis(epsilon, rep(0L, d), Gamma)

  # Quantile residuals
  rq <- as.numeric(stats::qnorm(pmin(pmax(
    get(paste0("maha_", copula))(mahalanobis, delta, d), EPS), 1 - EPS)))

  res <- switch(type,
                "quantile" = rq,
                "marginal" = rm,
                "epsilon" = epsilon)

  res

}

# Print
#' @rdname multireg-methods
#' @export
print.multireg <- function(x, ...)
{

  y <- if(is.null(x$y)){
    stats::model.response(stats::model.frame(x))
  } else{
    x$y
  }

  X <- stats::model.matrix(x)

  n <- x$nobs
  d <- x$d

  mu.link <- x$mu.link
  margins <- x$margins
  association <- x$association
  gamma <- x$gamma
  copula <- x$copula
  delta <- x$delta

  beta <- x$coefficients

  lambda_id <- apply(matrix(margins, ncol = 1), 1, function(x) !grepl("Log-", as.bcs(x)$name))
  nu_id <- apply(matrix(margins, ncol = 1), 1, function(x) as.bcs(x)$extrap)

  rownames(beta) <- colnames(X)
  colnames(beta) <- colnames(y)

  # Names
  copula_name <- paste(toupper(substr(copula, 1, 1)), substr(copula, 2, nchar(copula)), sep = "")
  gamma_name <- paste(toupper(substr(association, 1, 1)), substr(association, 2, nchar(association)), sep = "")

  if (requireNamespace("crayon", quietly = TRUE)) {

    cat(crayon::cyan(
      "\nMultivariate Box-Cox Regression with",
      if (is.null(delta)) copula_name else paste0(copula_name, "(", round(delta, 2), ")"),
      "Copula\n\n"
    ))

    cat(crayon::cyan("Call:\n\n"))
    print(x$call)
    if (x$optim_params$convergence > 0) {
      cat("\nmodel did not converge\n")
    } else {

      cat(crayon::cyan("\nRegression coefficients with", mu.link, "link function:\n\n"))
      print(round(beta, 4))

      tab <- cbind(sigma = x$fitted.values$sigma,
                   if(any(lambda_id)) lambda = x$fitted.values$lambda else NULL,
                   if(any(nu_id)) nu = x$fitted.values$nu else NULL)

      colnames(tab) <- c("sigma", if(any(lambda_id)) "lambda", if(any(nu_id)) "nu")

      cat(crayon::cyan("\nFurther marginal parameters:\n\n"))
      print(tab, na.print = "-")

      if (length(gamma) > 0) {

        cat(crayon::cyan("\n", gamma_name, " association matrix:\n\n", sep = ""))

        if (x$association == "non-associative"){
          Gamma <- diag(d)
        }else{
          Gamma <- get(tolower(x$association))(d)$Gamma(round(gamma, 4))
        }

        Gamma[upper.tri(Gamma)] <- diag(Gamma) <- NA
        colnames(Gamma) <- rownames(Gamma) <- colnames(y)
        print(Gamma, na.print = ".")
      }

      cat(
        "\n---",
        crayon::cyan("\nMargins:"), margins,
        crayon::cyan("\nlogLik:"), x$logLik, "|",
        crayon::cyan("AIC:"), stats::AIC(x), "|",
        crayon::cyan("BIC:"), stats::AIC(x, k = log(n)), "\n"
      )
    }

  } else {

    cat(
      "\nMultivariate Box-Cox Regression with",
      if (is.null(delta)) copula_name else paste0(copula_name, "(", round(delta, 2), ")"),
      "Copula\n\n"
    )

    cat("Call:\n\n")
    print(x$call)
    if (x$optim_params$convergence > 0) {
      cat("\nmodel did not converge\n")
    } else {

      cat("\nRegression coefficients with", mu.link, "link function:\n\n")
      print(round(beta, 4))

      tab <- cbind(sigma = x$fitted.values$sigma,
                   if(any(lambda_id)) lambda = x$fitted.values$lambda else NULL,
                   if(any(nu_id)) nu = x$fitted.values$nu else NULL)

      colnames(tab) <- c("sigma", if(any(lambda_id)) "lambda", if(any(nu_id)) "nu")

      cat("\nFurther marginal parameters:\n\n")
      print(tab, na.print = "-")

      if (length(gamma) > 0) {

        cat("\n", gamma_name, " association matrix:\n\n", sep = "")

        if (x$association == "non-associative"){
          Gamma <- diag(d)
        }else{
          Gamma <- get(tolower(x$association))(d)$Gamma(round(gamma, 4))
        }

        Gamma[upper.tri(Gamma)] <- diag(Gamma) <- NA
        colnames(Gamma) <- rownames(Gamma) <- colnames(y)
        print(Gamma, na.print = ".")
      }

      cat(
        "\n---",
        "\nMargins:", margins,
        "\nlogLik:", x$logLik, "|",
        "AIC:", stats::AIC(x), "|",
        "BIC:", stats::AIC(x, k = log(n)), "\n"
      )
    }

  }


  invisible(x)
}

# Summary
#' @rdname multireg-methods
#' @export
summary.multireg <- function(object, ...)
{

  y <- if(is.null(object$y)){
    stats::model.response(stats::model.frame(object))
  } else{
    object$y
  }

  X <- stats::model.matrix(object)

  n <- object$nobs
  d <- object$d
  b <- ncol(X)

  mu.link <- object$mu.link
  margins <- object$margins
  association <- object$association
  gamma <- object$gamma
  copula <- object$copula
  delta <- object$delta

  lambda_id <- apply(matrix(margins, ncol = 1), 1, function(x) !grepl("Log-", as.bcs(x)$name))
  nu_id <- apply(matrix(margins, ncol = 1), 1, function(x) as.bcs(x)$extrap)
  gamma_id <- length(gamma) > 0

  ## Summary for quantile residuals
  res <- stats::residuals(object)
  skewness <- mean((res - mean(res))^3) / (stats::sd(res)^3)
  kurtosis <- mean((res - mean(res))^4) / (stats::sd(res)^4)
  residuals_tab <- round(cbind(mean(res), stats::sd(res), skewness, kurtosis), 6)
  colnames(residuals_tab) <- c("Mean", "Std. dev.", "Skewness", "Kurtosis")
  rownames(residuals_tab) <- " "

  # Summary for mu
  beta_coef <- object$coefficients
  beta_se <- matrix(sqrt(diag(stats::vcov(object)[1:(b * d), 1:(b*d)])), ncol = d)
  beta_z <- beta_coef / beta_se
  beta_pvalue <- 2 * stats::pnorm(abs(beta_z), lower.tail = FALSE)

  rownames(beta_coef) <- rownames(beta_se) <- rownames(beta_z) <- rownames(beta_pvalue) <- colnames(X)
  colnames(beta_coef) <- colnames(beta_se) <- colnames(beta_z) <- colnames(beta_pvalue) <- colnames(y)

  # Summary for sigma
  sigma_coef <- sigma_se <- rep(NA, d)

  sigma_coef <- object$fitted.values$sigma
  sigma_se <- sqrt(diag(as.matrix(stats::vcov(object, "sigma"))))

  sigma_tab <- round(matrix(c(sigma_coef, sigma_se), byrow = TRUE, nrow = 2), 4)
  colnames(sigma_tab) <- colnames(y)
  rownames(sigma_tab) <- c(" ", "Std. Error")

  # Summary for lambda
  if (any(lambda_id)){

    lambda_coef <- lambda_se <- rep(NA, d)

    lambda_coef[lambda_id] <- object$fitted.values$lambda[lambda_id]
    lambda_se[lambda_id] <- sqrt(diag(as.matrix(stats::vcov(object, "lambda"))))

    lambda_tab <- round(matrix(c(lambda_coef, lambda_se), byrow = TRUE, nrow = 2), 4)
    colnames(lambda_tab) <- colnames(y)
    rownames(lambda_tab) <- c(" ", "Std. Error")

  } else {

    lambda_tab <- NULL

  }

  if (any(nu_id)) {

    nu_coef <- nu_se <- rep(NA, d)

    nu_coef[nu_id] <- object$fitted.values$nu[nu_id]
    nu_se[nu_id] <- sqrt(diag(as.matrix(stats::vcov(object, "nu"))))

    nu_tab <- round(matrix(c(nu_coef, nu_se), byrow = TRUE, nrow = 2), 4)
    colnames(nu_tab) <- colnames(y)
    rownames(nu_tab) <- c(" ", "Std. Error")

  } else {

    nu_tab <- NULL

  }

  # Summary for gamma
  gamma_tab <- NULL
  if (gamma_id) {

    if (object$association == "non-associative"){
      Gamma <- diag(d)
    }else{
      Gamma <- get(tolower(object$association))(d)$Gamma(round(gamma, 4))
    }

    Gamma[upper.tri(Gamma)] <- diag(Gamma) <- NA
    colnames(Gamma) <- rownames(Gamma) <- colnames(y)

  }

  out <- list(call = object$call,
              residuals = residuals_tab,
              mu.link = object$mu.link,
              beta_coef = beta_coef,
              beta_se = beta_se,
              beta_z = beta_z,
              beta_pvalue = beta_pvalue,
              sigma_tab = sigma_tab,
              lambda_tab = lambda_tab,
              nu_tab = nu_tab,
              Gamma = Gamma,
              margins = margins,
              association = association,
              copula = copula,
              delta = delta,
              logLik = object$logLik,
              n = n,
              d = d,
              y = y,
              lambda_id = lambda_id,
              nu_id = nu_id,
              gamma_id = gamma_id,
              AIC = stats::AIC(object),
              BIC = stats::AIC(object, k = log(n)))

  class(out) <- "summary.multireg"
  out
}

# Print summary
#' @rdname multireg-methods
#' @export
print.summary.multireg <- function(x, ...)
{

  y <- x$y
  n <- x$n
  d <- x$d
  copula <- x$copula
  delta <- x$delta
  association <- x$association
  margins <- x$margins
  mu.link <- x$mu.link

  lambda_id <- apply(matrix(margins, ncol = 1), 1, function(x) !grepl("Log-", as.bcs(x)$name))
  nu_id <- apply(matrix(margins, ncol = 1), 1, function(x) as.bcs(x)$extrap)

  copula_name <- paste(toupper(substr(copula, 1, 1)), substr(copula, 2, nchar(copula)), sep = "")
  gamma_name <- paste(toupper(substr(association, 1, 1)), substr(association, 2, nchar(association)), sep = "")

  if (requireNamespace("crayon", quietly = TRUE)) {

    cat(crayon::cyan(
      "\nMultivariate Box-Cox Regression with",
      if (is.null(delta)) copula_name else paste0(copula_name, "(", round(delta, 2), ")"),
      "Copula\n\n"
    ))

    cat(crayon::cyan("Call:\n\n"))
    print(x$call)

    cat(crayon::cyan("\nQuantile residuals:\n\n"))
    print(x$residuals)

    for (j in 1:d) {

      cat("\n", crayon::cyan(colnames(y)[j]), crayon::cyan(" ~ "),
          crayon::cyan(get(x$margins[j])()$name), crayon::cyan(" Distribution:\n"), sep = "")


      cat("\nRegression coefficients with", mu.link, "link function:\n")

      cmat <- round(cbind(Est = x$beta_coef[, j],
                    `Std. Error` = x$beta_se[, j],
                    `z value` = x$beta_z[, j],
                    `Pr(>|z|)` = x$beta_pvalue[, j]), 4)
      stats::printCoefmat(cmat)

      TAB <- rbind(
        x$sigma_tab[, j],
        if (lambda_id[j]) x$lambda_tab[, j],
        if (nu_id[j]) x$nu_tab[, j]
      )

      rownames(TAB) <- c("sigma", if (lambda_id[j]) "lambda", if (nu_id[j]) "nu")
      colnames(TAB) <- c("Est", "Std. Error")
      cat("\nFurther marginal parameters:\n")
      print(TAB)

    }

    if (x$gamma_id) {

      cat(crayon::cyan("\n", gamma_name, " association matrix:\n\n", sep = ""))
      print(x$Gamma, na.print = ".")

    }

    cat(
      "\n---",
      crayon::cyan("\nlogLik:"), x$logLik, "|",
      crayon::cyan("AIC:"), x$AIC, "|",
      crayon::cyan("BIC:"), x$BIC, "\n"
    )

  } else {


    cat(
      "\nMultivariate Box-Cox Regression with",
      if (is.null(delta)) copula_name else paste0(copula_name, "(", round(delta, 2), ")"),
      "Copula\n\n"
    )

    cat("Call:\n\n")
    print(x$call)

    cat("\nQuantile residuals:\n\n")
    print(x$residuals)

    for (j in 1:d) {

      cat("\n", colnames(y)[j], " ~ ",
          get(x$margins[j])()$name, " Distribution:\n", sep = "")


      cat("\nRegression coefficients with", mu.link, "link function:\n")

      cmat <- round(cbind(Est = x$beta_coef[, j],
                          `Std. Error` = x$beta_se[, j],
                          `z value` = x$beta_z[, j],
                          `Pr(>|z|)` = x$beta_pvalue[, j]), 4)
      stats::printCoefmat(cmat)

      TAB <- rbind(
        x$sigma_tab[, j],
        if (lambda_id[j]) x$lambda_tab[, j],
        if (nu_id[j]) x$nu_tab[, j]
      )

      rownames(TAB) <- c("sigma", if (lambda_id[j]) "lambda", if (nu_id[j]) "nu")
      colnames(TAB) <- c("Est", "Std. Error")
      cat("\nFurther marginal parameters:\n")
      print(TAB)

    }

    if (x$gamma_id) {

      cat("\n", gamma_name, " association matrix:\n\n", sep = "")
      print(x$Gamma, na.print = ".")

    }

    cat(
      "\n---",
      "\nlogLik:", x$logLik, "|",
      "AIC:", x$AIC, "|",
      "BIC:", x$BIC, "\n"
    )


    if (x$gamma_id) {

      cat("\n", gamma_name, " association matrix:\n\n", sep = "")
      print(x$Gamma, na.print = ".")

    }

    cat(
      "\n---",
      "\nlogLik:", x$logLik, "|",
      "AIC:", x$AIC, "|",
      "BIC:", x$BIC, "\n"
    )

  }


  invisible(x)
}



globalVariables(c("z", "upper", "lower", "density", "theo", "emp", "grid_x", "grid_y", "prob"))
# Plot
#' Visualization of the fit of the BCNSM distributions
#'
#' Five plots (selectable by \code{type}) are currently available:
#' a plot of residuals against fitted values, a plot of residuals against
#' the indexes, a Normal Q-Q plot, a barplot with comparisons of the
#' observed and fitted frequencies, and a plot of the sample autocorelations
#' of the residuals.
#'
#' @param x an object of class \code{multireg}.
#' @param type character; specifies which graphical should be produced. The available options
#'     are \code{"quantile"} (default), \code{"marginal"}, \code{"mahalanobis"}, and
#'     \code{"epsilon"}.
#' @param levels levels for contours plots. Used only when \code{type = "epsilon"}.
#' @param panel a vector of the form \code{c(nr, nc)} with the number of rows and columns, where
#'     the figures will be drawn in an \code{nr-}by\code{-nc} array on the device. Used
#'     only when \code{type = "marginal"}.
#' @param ... further arguments passed to or from other methods.
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @export
#'
plot.multireg <- function(x, type = c("quantile", "marginal", "mahalanobis", "epsilon"),
                          levels = c(1e-1, 1e-2, 1e-3, 1e-4), panel = NULL, ...) {

  type <- match.arg(type, c("quantile", "marginal", "mahalanobis", "epsilon"))

  y <- if(is.null(x$y)){
    stats::model.response(stats::model.frame(x))
  } else{
    x$y
  }

  n <- x$nobs
  d <- x$d

  copula <- x$copula
  delta <- x$delta
  Gamma <- get(x$association)(d)$Gamma(x$gamma)

  # Univariate marginal density of the NSM distribution
  if (copula == "slash") {
    dPSI <- function(x, log = FALSE) dslash(x, nu = delta, log = log)
  } else if (copula == "gaussian") {
    dPSI <- function(x, log = FALSE) stats::dnorm(x, log = log)
  } else if (copula == "t") {
    dPSI <- function(x, log = FALSE) stats::dt(x, delta, log = log)
  } else {
    dPSI <- function(x, log = FALSE) dhyp(x, nu = delta, log = log)
  }

  dmv <- get(paste0("dmv_", copula))

  # Information for the normal probability plot with a confidence region
  Pi <- (1:n - 0.5) / n
  zi <- stats::qnorm(Pi)

  if (is.null(panel)) panel <- c(ceiling(d / 4), 4)

  ### Quantile -------------------------------------------------------------------------------------
  if (type == "quantile") {

    rq <- sort(stats::residuals(x))
    mu_hat <- stats::median(rq)
    sigma_hat <- stats::IQR(rq) / 1.349

    se <- sigma_hat * sqrt(Pi * (1 - Pi) / n) / stats::dnorm(zi)
    rq_hat <- mu_hat + sigma_hat * zi

    positions <- data.frame(x = c(zi, rev(zi)),
                            y = c(rq_hat - stats::qnorm(1 - 0.05 / 2) * se,
                                  rev(rq_hat + stats::qnorm(1 - 0.05 / 2) * se)))

    ggplot2::ggplot(data.frame(zi = zi, rq = rq)) +
      ggplot2::geom_polygon(data = positions, ggplot2::aes(x = x, y = y), fill = "#cceff1") +
      ggplot2::geom_point(ggplot2::aes(x = zi, y = rq), pch = "+", cex = 4) +
      ggplot2::xlab("Normal quantiles") + ggplot2::ylab("Residual quantiles") +
      ggplot2::geom_qq_line(ggplot2::aes(sample = rq), col = "#56B1F7", lwd = 1)

    ### Marginal -----------------------------------------------------------------------------------
  } else if (type == "marginal") {

    rm <- stats::residuals(x, "marginal")

    res <- data.frame(emp = c(apply(rm, 2, sort)),
                  theo = rep(stats::qnorm(stats::ppoints(n)), d),
                  marg = factor(rep(colnames(y), each = n), levels = colnames(y)))

    lower_f <- function(res, z) {
      robust_sd <- stats::IQR(res) / 1.349
      robust_line <- stats::median(res) + z * robust_sd
      robust_se <-  robust_sd * sqrt(stats::pnorm(z) * (1 - stats::pnorm(z)) / n) / stats::dnorm(z)
      robust_line - stats::qnorm(1 - 0.05/2) * robust_se
    }

    upper_f <- function(res, z) {
      robust_sd <- stats::IQR(res) / 1.349
      robust_line <- stats::median(res) + z * robust_sd
      robust_se <-  robust_sd * sqrt(stats::pnorm(z) * (1 - stats::pnorm(z)) / n) / stats::dnorm(z)
      robust_line + stats::qnorm(1 - 0.05/2) * robust_se
    }

    band <- data.frame(
      lower = c(apply(rm, 2, lower_f, z = seq(min(res$theo), max(res$theo), length.out = 500))),
      upper = c(apply(rm, 2, upper_f, z = seq(min(res$theo), max(res$theo), length.out = 500))),
      marg = factor(rep(colnames(y), each = 500), levels = colnames(y)),
      z = rep(seq(min(res$theo), max(res$theo), length.out = 500), d)
    )

    ggplot2::ggplot(res) +
      ggplot2::geom_ribbon(ggplot2::aes(x = z, ymax = upper, ymin = lower),  data = band, fill = "#cceff1",  color = "#cceff1",
                  show.legend = FALSE) +
      ggplot2::geom_point(ggplot2::aes(x = theo, y = emp), pch = "+", cex = 2) +
      ggplot2::geom_qq_line(ggplot2::aes(sample = emp), col = "#56B1F7", lwd = 1) +
      ggplot2::facet_wrap(~ marg, nrow = panel[1], ncol = panel[2]) +
      ggplot2::labs(x = "Theoretical quantiles",
           y = "Sample quantiles")

    ### Mahalanobis ----------------------------------------------------------------------------------
  } else if (type == "mahalanobis") {

    epsilon <- stats::residuals(x, type = "epsilon")
    mahalanobis <- mahalanobis(epsilon, rep(0L, d), Gamma)

    ggplot2::ggplot() +
      ggplot2::geom_segment(ggplot2::aes(x = 1:n, y = rep(0, n),
                       xend = 1:n, yend = mahalanobis)) +
      ggplot2::labs(x = "Index observations", y = "Transformed Mahalanobis distances")

    ### Epsilon --------------------------------------------------------------------------------------
  } else if (type == "epsilon") {

    epsilon <- stats::residuals(x, "epsilon")

    diag_func <- function(data, mapping, ...){

      x <- GGally::eval_data_col(data, mapping$x)

      ggplot2::ggplot(data, mapping) +
        ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                                colour = 1, fill = "white",
                                bins = ceiling(1 + 3.33 * log(n))) +
        ggplot2::geom_function(fun = dPSI, col = "#00AFBB")


    }

    ### Upper plots
    upper_func <- function(data, mapping, ...){

      x <- GGally::eval_data_col(data, mapping$x)
      y <- GGally::eval_data_col(data, mapping$y)

      i <- which(colSums(epsilon - x) == 0)
      j <- which(colSums(epsilon - y) == 0)

      corr <- Gamma[i, j]

      colFn <- grDevices::colorRampPalette(c("brown1", "white", "dodgerblue"), interpolate ='spline')
      fill <- colFn(1000)[findInterval(corr, seq(-1, 1, length = 1000))]

      GGally::ggally_text(
        label = as.character(round(corr, 4)),
        mapping = ggplot2::aes(),
        xP = 0.5,
        yP = 0.5,
        cex = 2.5,
        color = 'black', ...) +
        ggplot2::theme_void() +
        ggplot2::theme(panel.background = ggplot2::element_rect(fill = fill))
    }

    ### Lower plots
    lower_func <- function(data, mapping, ...){

      x <- GGally::eval_data_col(data, mapping$x)
      y <- GGally::eval_data_col(data, mapping$y)

      i <- which(colSums(epsilon - x) == 0)
      j <- which(colSums(epsilon - y) == 0)

      Gamma_aux <- matrix(c(1, Gamma[i, j], Gamma[i, j], 1), 2, 2)


      grid <- expand.grid(seq(min(x) - 10, max(x) + 10, length.out = 200),
                          seq(min(y) - 10, max(y) + 10, length.out = 200))

      data_aux <- data.frame(grid_x = grid[, 1],
                             grid_y = grid[, 2],
                             prob = dmv(cbind(grid[, 1], grid[, 2]),
                                        mu = rep(0, 2L),
                                        Sigma = Gamma_aux,
                                        delta = delta))


      ggplot2::ggplot() +
        ggplot2::geom_point(data = data, mapping = ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_contour(data = data_aux, mapping =  ggplot2::aes(x = grid_x, y = grid_y, z = prob),
                              col = "#00AFBB", breaks = levels) +
        ggplot2::labs(x = "", y = "")


    }

    colFn <- grDevices::colorRampPalette(c("brown1", "white", "dodgerblue"), interpolate ='spline')
    lab <- ggplot2::ggplot(data.frame(x = stats::runif(200, -1, 1), y = stats::runif(200, -1, 1),
                                      z = stats::runif(200, -1, 1)), ggplot2::aes(x, y, colour = z)) +
      ggplot2::geom_point() +
      ggplot2::scale_colour_gradient2("Fitted association \nparameter",
                                      low = colFn(200)[1], high = colFn(200)[200]) +
      ggplot2::theme(legend.title.align = 0.5,
                     legend.position = "top",
                     legend.key.height = ggplot2::unit(0.3, 'cm'),
                     legend.key.width = ggplot2::unit(1.5, 'cm'))

    GGally::ggpairs(as.data.frame(epsilon),
                    upper = list(continuous = GGally::wrap(upper_func)),
                    lower = list(continuous = GGally::wrap(lower_func)),
                    diag = list(continuous = GGally::wrap(diag_func)),
                    legend = GGally::grab_legend(lab),
                    progress = FALSE) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::theme(legend.position = "top")

  }


}
