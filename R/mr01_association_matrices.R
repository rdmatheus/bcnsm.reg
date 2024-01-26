#' @name association
#'
#' @title Association Matrices for the BCNSM Marginal Regression Models
#'
#' @description Current available association matrices in the \code{bcnsm.reg}
#'     package for the class of the BCNSM marginal regression models.
#'
#' @param p,q Order of the autoregressive and the moving average component,
#'     respectively, passed as arguments to \code{arma()}.
#' @param id 	Subject id for longitudinal/clustered data. This is a vector of
#'     the same length of the number of observations. Please note that data
#'     must be sorted in way that observations from the same cluster are
#'     contiguous.
#' @param type  A character string specifying the correlation structure among
#'     groups for longitudinal/clustered data. At the
#'     moment, the following are implemented:
#'     \tabular{ll}{
#'     \code{nonassociative}  \tab Non-associative. \cr
#'     \code{ar1}  \tab Autoregressive of order 1. \cr
#'     \code{ma1}  \tab Moving average of order 1. \cr
#'     \code{uniform1}  \tab Uniform or exchangeable. \cr
#'     \code{unstructured}  \tab Unstructured \cr
#'     }
#' @param D Matrix with values of the distances between pairs of data
#'     locations for spatial data.
#' @param smt Value of the shape parameter of the Matérn correlation class.
#'     The default \code{smt = 0.5} corresponds to an exponential correlation
#'     model.
#'
#' @details
#'
#' The functions related to the association structures are a direct
#'     adaptation of the functions of the \code{\link[gcmr]{gcmr}} package.
#'     The documentation of the original functions of the package can be seen
#'     at \code{\link[gcmr]{cormat.gcmr}} documentation. The available
#'     association matrices are:
#'  \tabular{ll}{
#'  \bold{Function}  \tab \bold{Dependence structure}\cr
#'  \code{nonassociative}  \tab Non-associative structure. \cr
#'  \code{arma}  \tab ARMA(p, q). \cr
#'  \code{cluster}  \tab Longitudinal/clustered data. \cr
#'  \code{matern}  \tab Matérn spatial correlation. \cr
#'  }
#'
#'  @return
#'
#'  Each type of structure will require different arguments for the association
#'  matrix to be completed. However, all functions return a list of the
#'  following components:
#'  \itemize{
#'  \item{\code{npar}:}{ Number of parameters associated with the correlation
#'      structure.}
#'  \item{\code{start}:}{ A function \code{function(y)} that returns the
#'      initial value for use in optimization algorithms. Its argument is the
#'      vector/matrix of observations \code{y}.}
#'  \item{\code{Gamma}:}{ A function \code{function(alpha, d)} which returns the
#'      association matrix itself. Its arguments are the parameter vector
#'      associated with the correlation matrix structure and the dimension
#'      of the correlation matrix.}
#'  \item{\code{name}:}{ Name of the specified correlation structure.}
#'  }
#'
#' @references Guido Masarotto, & Cristiano Varin (2017). Gaussian Copula
#'     Regression in R. Journal of Statistical Software, 77(8), 1--26.
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @export

## ARMA(p,q) working correlation for time-series ---------------------------------------------------
#' @rdname association
#' @export
arma <- function(p = 0L, q = 0L) {

  ar <- if(p) 1:p else NULL
  ma <- if(q) 1:q + p else NULL

  ans <- list()

  ans$npar <- p + q

  ans$start <- function(y = NULL) {

    n <- length(y)

    if (p | q){
      gamma <- rep(0, p + q)
      names(gamma) <- c(if (p) paste0("ar", 1:p) else NULL, if (q) paste0("ma", 1:q) else NULL)
    } else {
      gamma <- NULL
    }

    gamma
  }

  ans$Gamma <- function(gamma, n){

    if (p | q){

      if((p && any(Mod(polyroot(c(1,-gamma[ar]))) < 1.01) ) || (q && any(Mod(polyroot(c(1, gamma[ma]))) < 1.01)))
        return(NULL)

      gamma <- stats::ARMAacf(gamma[ar], gamma[ma], n - 1)
      r <- seq(1, n)
      Matrix::as.matrix(outer(r, r, function(i, j) gamma[1 + abs(i - j)]))

    } else {

      diag(n)

    }


  }

  ans$name <- paste0("ARMA(", p, ", ", q,")")

  ans
}

## Matern working correlation for spatial data -----------------------------------------------------

## D is a distance matrix, smt is the smoothing parameter
#' @rdname association
#' @export
matern <- function(D, smt = 0.5) {
  ans <- list()
  ans$npar <- 1
  ans$start <- function(y = NULL) {
    gamma <- stats::median(D)
    names(gamma) <- c("gamma")
    attr(gamma,"lower") <- sqrt(.Machine$double.eps)
    gamma
  }

  ans$Gamma <- function(gamma, n){
    S <- try(Matrix::forceSymmetric(.matern(D, gamma, smt)), silent = TRUE)
    if(inherits(S, "try-error") | gamma < 0) NULL else S
  }

  ans$name <- "Matern spatial"
  ans
}

# Adapted from 'geoR' package
.matern <- function (u, phi, kappa)
{
  if (is.vector(u))
    u <- matrix(u)

  dimnames(u) <- list(NULL, NULL)

  uphi <- u/phi

  out <- Matrix::Matrix(NaN, dim(u)[1], dim(u)[2])

  id1 <- which(u > 0 & phi > 0 & kappa > 0, arr.ind = TRUE)
  id2 <- which(u <= 0 & phi > 0 & kappa > 0, arr.ind = TRUE)

  out[id1] <- exp(-(kappa - 1) * log(2) - lgamma(kappa) + kappa * log(uphi[id1]) + log(besselK(uphi[id1], kappa)))
  out[id2] <- 1
  out[which(u > 600 * phi, arr.ind = TRUE)] <- 0

  out
}

### Longitudinal/Clustered data working correlation ----------------------------
# It assumes that it is not possible that all the observations inside a cluster
# can be missing
#' @rdname association
#' @export
cluster <- function(id, type = c("nonassociative", "ar1", "ma1",
                                 "uniform", "unstructured")) {

  type <- match.arg(type, c("nonassociative", "ar1", "ma1",
                            "uniform", "unstructured"))

  if(!length(rle(id)$values) == length(unique(id)))
    stop("data must be sorted in way that observations from the same cluster are contiguous")

  ng <- 1:length(unique(id))
  if (!(length(ng) > 1)) stop("only one strata")

  if (type == "nonassociative") {
    ans <- nonassociative(length(id))
    ans$id <- id
    return(ans)
  }

  ans <- list(type = type, id = id)
  ans$npar <- if (type != "unstructured") 1 else choose(max(table(id)), 2)
  data <- data.frame(id = id)
  fn <- switch(type,
               "ar1" = function(g) nlme::corAR1(g, form = ~ 1 | id),
               "ma1" = function(g) nlme::corARMA(g, form = ~ 1 |id, p = 0, q = 1),
               "uniform" = function(g) nlme::corCompSymm(g, form = ~ 1 | id),
               "unstructured" = function(g) nlme::corSymm(g, form = ~ 1 | id))
  ans$start <- function(y = NULL) {
    np <-  if(type != "unstructured") 1 else choose(max(table(id)), 2)
    alpha <- rep(0, np)
    names(alpha) <- switch(type, "ar1" = "ar1", "ma1" = "ma1",
                           "uniform" = "alpha",
                           "unstructured" = paste("alpha", 1:ans$npar, sep = ""))
    eps <- sqrt(.Machine$double.eps)
    attr(alpha, "lower") <- rep(-1 + eps, np)
    attr(alpha, "upper") <- rep(1 - eps, np)
    alpha
  }
  ans$Gamma <- function(alpha, n) {
    q <- try(nlme::corMatrix(nlme::Initialize(fn(alpha), data = data)),
             silent = TRUE)
    if (inherits(q, "try-error")) return(NULL)
    g <- split(rep(TRUE, n), id)
    q <- try(lapply(ng, function(i) q[[i]][g[[i]],g[[i]]]), silent = TRUE)
    if (inherits(q, "try-error") ) NULL else as.matrix(Matrix::bdiag(q))
  }
  ans$name <- "Clustered"
  ans
}

