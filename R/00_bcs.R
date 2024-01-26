#' @name bcs
#' @title Methods for \code{"bcs"} objects
#' @param x,object an object of class \code{"bcs"}.
#' @param ... further arguments passed to or from other methods.
#'
#' @details The class of the Box-Cox symmetric (BCS) distributions was introduced by
#'     Ferrari and Fumes (2017). The class includes the Box-Cox t (Rigby and Stasinopoulos, 2006),
#'     Box-Cox normal (or Box-Cox Cole-Green; Cole and Green, 1992), Box-Cox power exponential
#'     (Rigby and Stasinopoulos, 2004) distributions, and the class of the
#'     log-symmetric distributions (Vanegas and Paula, 2016) as special cases.
#'     The current available BCS distributions can be seen below.
#'   \tabular{llc}{
#'  \bold{Distribution}  \tab \bold{Abbreviation} \tab \bold{N. of parameters}\cr
#'  Box-Cox Cauchy  \tab \code{"bcca"}      \tab  3  \cr
#'  Box-Cox Hyperbolic  \tab \code{"bchp"}      \tab  4  \cr
#'  Box-Cox Laplace  \tab \code{"bcla"}      \tab  3  \cr
#'  Box-Cox Type I Logistic  \tab \code{"bcloi"}      \tab  3  \cr
#'  Box-Cox Type II Logistic  \tab \code{"bcloii"}      \tab  3  \cr
#'  Box-Cox Normal  \tab \code{"bcno"}      \tab  3  \cr
#'  Box-Cox Power exponential  \tab \code{"bcpe"}      \tab  4  \cr
#'  Box-Cox Slash  \tab \code{"bcsl"}      \tab  4  \cr
#'  Box-Cox t  \tab \code{"bct"}      \tab  4  \cr
#'  }
#'
#'  Log-symmetric special cases are also available:
#'   \tabular{llc}{
#'  \bold{Distribution}  \tab \bold{Abbreviation} \tab \bold{N. of parameters}\cr
#'  Log-Cauchy  \tab \code{"lca"}      \tab  2  \cr
#' Log-Hyperbolic  \tab \code{"lhp"}      \tab  3  \cr
#'  Log-Laplace  \tab \code{"lla"}      \tab  2  \cr
#'  Log-Type I Logistic  \tab \code{"lloi"}      \tab  2  \cr
#'  Log-Type II Logistic  \tab \code{"lloii"}      \tab  2  \cr
#'  Log-Normal  \tab \code{"lno"}      \tab  2  \cr
#'  Log-Power exponential  \tab \code{"lpe"}      \tab  3  \cr
#' Log-Slash  \tab \code{"lsl"}      \tab  3  \cr
#'  Log-t  \tab \code{"lt"}      \tab  3  \cr
#'  }
#'
#'
#' @references
#'  Cole, T., and Green, P.J. (1992). Smoothing reference centile curves: the LMS
#'      method and penalized likelihood. \emph{Statistics in medicine}, 11, 1305-1319.
#'
#'  Rigby, R. A., and Stasinopoulos, D. M. (2004). Smooth centile curves for skew
#'      and kurtotic data modelled using the Box-Cox power exponential
#'      distribution. \emph{Statistics in medicine}, 23, 3053-3076.
#'
#'  Rigby, R. A., and Stasinopoulos, D. M. (2006). Using the Box-Cox t
#'      distribution in GAMLSS to model skewness and kurtosis. \emph{Statistical Modelling}, 6, 209-229.
#'
#'  Vanegas, L. H., and Paula, G. A. (2016). Log-symmetric distributions:
#'      statistical properties and parameter estimation. \emph{Brazilian Journal of Probability and Statistics}, 30, 196-220.
#'
#'  Ferrari, S. L., and Fumes, G. (2017). Box-Cox symmetric distributions and
#'      applications to nutritional data. \emph{AStA Advances in Statistical Analysis}, 101, 321-344.
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @examples
#' bcs("bcno")
#' bcs("lno")
#'
#' bcs("bct")
#' bcs("lt")
#'
NULL

# Class definition
#' @rdname bcs
#' @export
bcs <- function(object, ...) UseMethod("bcs")

# Default
#' @rdname bcs
#' @export
bcs.default <- function(object, ...) {
  cl <- data.class(object)[1]

  return(switch(cl,
                bcs = object,
                "function" = bcs(object()),
                character = bcs(get(object)),
                name = bcs(eval(object)),
                call = bcs(eval(object)),
                NULL = bcno(),
                stop("The object argument is invalid")
  ))
}

# Transformation to bcs class
#' @rdname bcs
#' @export
as.bcs <- function(object) {
  if (inherits(object, "bcs")) object else bcs(object)
}

# Print method
#' @rdname bcs
#' @export
print.bcs <- function(x, ...) {

  if (requireNamespace("crayon", quietly = TRUE)) {

    cat(
      crayon::cyan("The class of the Box-Cox symmetric distributions\n"),
      crayon::cyan("\nName:"), x$name,
      crayon::cyan("\nAbbreviation:"), x$abb,
      crayon::cyan("\nNumber of parameters:"), x$npar, "\n"
    )

  } else {

    cat(
      "The class of the Box-Cox symmetric distributions\n",
      "\nName:", x$name,
      "\nAbbreviation:", x$abb,
      "\nNumber of parameters:", x$npar, "\n"
    )

  }



}







