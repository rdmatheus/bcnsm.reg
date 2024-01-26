#' BCNSM Marginal Regression Models for Positive Continuous Data
#'
#' Fit a BCNSM marginal regression model to correlated positive data with Box-Cox symmetric (BCS)
#'  distributions.
#'
#'
#' @param formula a simbolic description of the model, of type
#'     \code{y ~ x} for covariates in the mean submodel only or \code{y ~ x | z}
#'     to enter covariates in the variance submodel. See details below.
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
#' @param margin a character which specifies the marginal BCS distribution.
#'     A table with the current available BCS distributions can be seen in \code{\link{bcs}}.
#' @param mu.link,sigma.link character; specifies the link function for the \code{mu} and
#'     \code{sigma} submodels, respectively. The links \code{"log"} (default) and
#'     \code{"identity"} are currently available.
#' @param association a function which returns an object of class \code{"\link{association}"},
#'     that describes the dependence structure of the data. See \code{\link{association}}.
#' @param copula character; informs which normal scale mixture distribution
#'     should be used to generate the copula. Currently,
#'     the copulas available are: Gaussian (\code{"gaussian"}),
#'     Student's t (\code{"t"}), slash (\code{"slash"}), and hyperbolic (\code{"hyp"}).
#' @param delta possible extra parameter associated with the copula. For example, the degrees of
#'     freedom of the \code{t} copula.
#' @param y,x logical; if \code{TRUE} the corresponding components of the fit, response and model matrix, are returned.
#' @param control  a list of optimization control arguments specified via \code{control_fit}  (under development).
#' @param ... further optimization control arguments passed to \code{control_fit} (under development).
#'
#' @details
#'
#' The \code{marginreg} function implements a broad class of marginal regression models for
#'     analyzing correlated positive data with Box-Cox symmetric marginal distributions.
#'     The dependence is described by a normal scale mixture copula, offering a fully-parametric
#'     alternative to the classical generalized estimating equation models. The proposed class
#'     provides a flexible modeling framework for positive data with different levels of skewness
#'     and tail-heaviness, allowing for different dependence structures, including independence,
#'     time series, longitudinal, clustered, or spatially correlated data through the use of
#'     the copula.
#'
#' The normal scale mixture copula entirely describes the dependence structure of the
#'     models. Because copulas are invariant under strictly increasing transformations,
#'     copula-based measures of association of the scale mixtures of normal distributions are the
#'     same for the BCNSM distributions. This statement is valid for measures of concordance, such
#'     as Kendall's tau or the tail dependence coefficients.
#'
#'     Let \code{gamma_ij} be the element of the ith row and jth column of the association matrix,
#'     say an n x n correlation matrix \code{Gamma}, where \code{gamma_ii = 1}, for \code{i, j = 1, ..., n}. If
#'     \code{Y = (Y1, ..., Yn)}' is an n-dimensional random vector with a BCNSM distribution and with
#'     association matrix \code{Gamma}. Thus, Kendall's tau of \code{Yi} and \code{Yj} is given by
#'
#'     \code{tau(Yi, Yj) = 2 * asin(gamma_ij) / pi},
#'
#'     for \code{i, j = 1, ..., n} and it is invariant under the copula. Since Kendall’s tau is a
#'     rank correlation coefficient, it measures any monotonous association between two random variables.
#'     Thus, we have an increasing monotonous association between \code{Yi} and \code{Yj} when \code{gamma_ij}
#'     is close to 1, and a decreasing monotonous association when \code{gamma_ij} is close to -1.
#'     Furthermore, \code{gamma_ij = 0} indicates no monotonous association. These facts show that
#'     the association matrix \code{Gamma} plays an essential role in specifying the dependence
#'     structure of \code{Y}. It can be parameterized in terms of a real-valued vector \code{gamma}
#'     in an interpretative way depending on the type of association that is of interest to describe.
#'     Currently, the \code{bcnsm.reg} package provides the non-associative, ARMA(p, q), cluster, and
#'     matern association structures.
#'
#'     The basic formula is of type \code{y ~ x1 + x2 + ... + xp} which
#'     specifies the model for the \code{mu} submodel only with \code{p} explanatory
#'     variables. Following the syntax of the \code{betareg} package
#'     (Cribari-Neto and Zeileis, 2010), the model for the \code{sigma} submodel, say in terms of
#'     z1, z2, ..., zk, is specified as \code{y ~ x1 + x2 + ... + xp | z1 + z2 +
#'     ... +zk} using functionalities inherited from package \code{Formula}
#'     (Zeileis and Croissant, 2010).
#'
#' @return The \code{marginreg} function returns an object of class "\code{marginreg}",
#'  which consists of a list with the following components:
#' \describe{
#'   \item{coefficients}{a list containing the elements \code{"beta"} and \code{"kappa"}.
#'       which consists of the estimates of the regression coefficients associated with the
#'       \code{mu} and the \code{sigma} submodels, respectively.}
#'   \item{fitted.values}{a list with the fitted values for \code{mu}, \code{sigma},
#'       \code{lambda}, and \code{nu}.}
#'   \item{margin}{a character with the marginal BCS distribution of the fit.}
#'   \item{links}{a list with elements \code{"mu.link"} and \code{"sigma.link"} with the specified link
#'       functions for the \code{mu} and \code{sigma} submodels, respectively.}
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
#'   \item{formula}{the formula used to specify the model in \code{marginreg}.}
#'   \item{terms}{the terms for the "\code{mu}" submodel.}
#'  }
#'
#' @export
#'
#' @references
#'
#' Cribari-Neto F, Zeileis A (2010). Beta regression in R. \emph{Journal
#'  of Statistical Software}, \bold{34}, 1–24
#'
#'  Ferrari, S. L., and Fumes, G. (2017). Box-Cox symmetric distributions and
#'      applications to nutritional data. \emph{AStA Advances in Statistical Analysis}, 101, 321-344.
#'
marginreg <- function(formula, data, subset, na.action, margin = "bcno",
                      mu.link = "log", sigma.link = "log",
                      association = nonassociative(),
                      copula = "gaussian", delta = NULL,
                      control = control_fit(...), y = FALSE, x = FALSE, ...)
{

  ret.y <- y
  ret.x <- x

  ## Call
  cl <- match.call()
  if (missing(data))  data <- environment(formula)

  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE

  ## Formula
  oformula <- stats::as.formula(formula)
  formula <- Formula::as.Formula(formula)

  if (length(formula)[2L] < 2L) {

    formula <- Formula::as.Formula(formula(formula), ~ 1)

  }else {

    if (length(formula)[2L] > 2L) {
      formula <- Formula::Formula(formula(formula, rhs = 1:2))
      warning("formula must not have more than two RHS parts")
    }

  }

  mf$formula <- formula

  ## Evaluate model.frame
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  ## Extract response and model matrices
  y <- stats::model.response(mf)
  X <- stats::model.matrix(formula, mf, rhs = 1L)
  Z <- stats::model.matrix(formula, mf, rhs = 2L)

  ## Some conditions
  if (length(y) < 1)
    stop("empty model\n")
  if (any(y < 0))
    stop("invalid dependent variable, all observations must be greater than zero\n")

  n <- length(y)

  if (is.null(association)) association <- nonassociative(n)

  ## Fit
  out <- mr_mle(y, X, Z, association, mu.link, sigma.link, copula, delta, margin, control)

  ## Further model information
  out$call <- cl
  out$formula <- formula
  out$terms <- list(mu = stats::terms(formula, data = data, rhs = 1L),
                    sigma = stats::terms(formula, data = data, rhs = 2L),
                    full = stats::terms(formula, data = data))

  if(ret.y) out$y <- y
  if(ret.x) out$x <- list(mu = X, sigma = Z)

  class(out) <- "marginreg"
  out
}
