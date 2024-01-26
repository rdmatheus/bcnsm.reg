#' BCNSM Regression for Multivariate Positive Data
#'
#' @description Fit of the BCNSM multivariate regression model via maximum
#' likelihood.
#'
#' @param formula a symbolic description of the model to be fitted. For example,
#'     \code{formula = cbind(y1, y2, y3) ~ x1 + x2} fits a BCNSM regression for the response
#'     variables y1, y2, y3 as a function of the explanatory variables x1 and x2.
#' @param data an optional data frame containing the variables in the formula. By default the
#'     variables are taken from environment(formula).
#' @param subset an optional vector specifying a subset of observations to be used in the fitting
#'     process. (See additional details about how this argument interacts with data-dependent bases
#'     in the ‘Details’ section of the \code{\link[stats]{model.frame}} documentation.)
#' @param na.action a function which indicates what should happen when the data contain \code{NAs}.
#'     The default is set by the \code{na.action} setting of \link[base]{options}, and is
#'     \code{\link[stats]{na.fail}} if that is unset. The ‘factory-fresh’ default is
#'     \code{\link[stats]{na.omit}}. Another possible value is \code{NULL}, no action.
#'     Value \code{\link[stats]{na.exclude}} can be useful.
#' @param margins a character or a character vector; specifies the marginal BCS distributions. If all
#'     BCS margins are the same, it is sufficient to enter only one character. A table with the
#'     current available BCS distributions can be seen in \code{\link{bcs}}.
#' @param mu.link character; specifies the link function for the \code{mu}
#'     submodel. The links \code{"log"} (default) and \code{"identity"} are currently available.
#' @param association one of \code{"unstructured"} (default), \code{"uniform"}, or
#'     \code{"nonassociative"}, which specify the association matrix of the model.
#' @param copula character; informs which normal scale mixture distribution
#'     should be used to generate the NSM copula. Currently,
#'     the copulas available are: Gaussian (\code{"gaussian"}),
#'     Student's t (\code{"t"}), slash (\code{"slash"}), and hyperbolic (\code{"hyp"}).
#' @param delta possible extra parameter associated with the copula. For example, the degrees of
#'     freedom of the \code{t} copula.
#' @param y,x logical; if \code{TRUE} the corresponding components of the fit, response and model matrix, are returned.
#' @param control  a list of optimization control arguments specified via \code{control_fit}  (under development).
#' @param ... further optimization control arguments passed to \code{control_fit} (under development).
#'
#' @return The \code{bcnsmgreg} function returns an object of class "\code{bcnsmreg}",
#'  which consists of a list with the following components:
#' \describe{
#'   \item{coefficients}{the matrix of regression coefficients associated with the \code{mu} submodel.}
#'   \item{fitted.values}{a list with the fitted values for \code{mu}, \code{sigma},
#'       \code{lambda}, and \code{nu}.}
#'   \item{margins}{a character vector with the marginal BCS distributions of the fit.}
#'   \item{mu.link}{specified link function for the \code{mu} submodel.}
#'   \item{copula, delta}{\code{"copula"} is a character which informs which normal scale mixture distribution
#'       was used to generate the NSM copula and \code{"delta"} is the possible extra parameter associated with
#'       the copula.}
#'   \item{gamma}{the estimated parameters of the association matrix, if any.}
#'   \item{association}{structure of the association matrix. It can be one of \code{"non-associative"},
#'       \code{"unstructured"}, or \code{"uniform"}.}
#'   \item{logLik}{log-likelihood of the fitted model.}
#'   \item{vcov}{asymptotic covariance matrix of the maximum likelihood estimator of the model parameters vector.
#'       Specifically, the inverse of the observed information matrix, obtained via numeric Hessian matrix.}
#'   \item{y}{the response vector (if \code{y = TRUE}).}
#'   \item{x}{the model matrix from the "\code{mu}" submodel, (if \code{x = TRUE}).}
#'   \item{optim_params}{control optimization parameters used by \code{\link{control_fit}}.}
#'   \item{nobs,d}{the number of observations in the sample and the dimension of the response variable, respectively.}
#'   \item{call}{ the function call.}
#'   \item{formula}{the formula used to specify the model in \code{bcnsmreg}.}
#'   \item{terms}{the terms for the "\code{mu}" submodel.}
#'  }
#'
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @examples
#' \dontrun{
#'
#' # Winscosin Breast Cancer Dataset
#' ?wdbc
#'
#' # Training set index
#' set.seed(123)
#' id <- sample(1:nrow(wdbc), 0.7 * nrow(wdbc))
#'
#' # Reference model
#' fit0 <- bcnsmreg(cbind(Texture, Area, Smoothness, Compactness, Concavity) ~ Diagnosis,
#'                  data = wdbc, subset = id)
#' fit0
#'
#' ## Marginal quantile residuals of the reference model
#' plot(fit0, "marginal", panel = c(2, 3))
#'
#' # Improve the fit on margins
#' fit_gaussian <- bcnsmreg(cbind(Texture, Area, Smoothness, Compactness, Concavity) ~ Diagnosis,
#'                          data = wdbc, subset = id, margins = c("lt", "lt", "lno", "lpe", "bct"))
#'
#' ## Marginal quantile residuals of the improved model
#' plot(fit_gaussian, "marginal", panel = c(2, 3))
#'
#'
#' ## Summary
#' summary(fit_gaussian)
#'
#' ## Overall quantile residuals of the final model
#' plot(fit_gaussian)
#'
#' ## The epsilon's transformations
#' plot(fit_gaussian, "epsilon")
#' }
#'
#'
bcnsmreg <- function(formula, data, subset, na.action,
                     margins = "bcno",
                     mu.link = "log",
                     association = c("unstructured", "uniform", "nonassociative"),
                     copula = c("gaussian", "t", "slash", "hyp"),
                     delta = NULL,
                     y = FALSE,
                     x = FALSE,
                     control = control_fit(...), ...)
{

  ret.y <- y
  ret.x <- x

  ## Call
  cl <- match.call()
  if (missing(formula)) stop("a formula argument is required")
  if (missing(data))  data <- environment(formula)

  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE

  ## Formula
  oformula <- stats::as.formula(formula)
  formula <- Formula::as.Formula(formula)

  if (length(formula)[2L] > 1L) {
    warning("formula must not have more than one RHS parts")
  }

  formula <- Formula::Formula(formula(formula, rhs = 1))

  mf$formula <- formula

  ## Evaluate model.frame
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  ## Extract terms, model matrix, response
  mt <- stats::terms(formula, data = data)
  mtX <- stats::terms(formula, data = data, rhs = 1L)

  y <- stats::model.response(mf, "numeric")
  n <- nrow(y)
  d <- ncol(y)

  X <- stats::model.matrix(formula, mf, rhs = 1L)
  b <- ncol(X)

  if (n < 1)
    stop("empty model")

  copula <- match.arg(copula, c("gaussian", "t", "slash", "hyp"))
  association <- match.arg(association, c("unstructured", "uniform", "nonassociative"))
  if(missing(delta)) delta <- NULL

  out <- bcnsm_mle(y = y, X = X, association = association,
                   copula = copula, delta = delta, margins = margins,
                   mu.link = mu.link, control = control)

  out$terms <- mtX

  if(ret.y) out$y <- y
  if(ret.x) out$x <- X

  ## Further model information
  out$call <- cl
  out$formula <- formula

  class(out) <- "multireg"
  out

}
