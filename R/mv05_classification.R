#' BCNSM Classification
#'
#' @param object an object of class \code{"multireg"}.
#' @param test matrix or data frame of test set cases.
#' @param cl factor of true classifications of training set.
#'
#' @return Factor of classifications of test set.
#'
#' @export
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
#' # BCNSM regression
#' fit <- bcnsmreg(cbind(Texture, Area, Smoothness, Compactness, Concavity) ~ Diagnosis,
#'                 data = wdbc, subset = id, margins = c("lt", "lt", "lno", "lpe", "bct"))
#'
#' ## Marginal quantile residuals
#' plot(fit_gaussian, "marginal", panel = c(2, 3))
#'
#' ## Summary
#' summary(fit_gaussian)
#'
#' ## Overall quantile residuals of the final model
#' plot(fit_gaussian)
#'
#' ## The epsilon's transformations
#' plot(fit_gaussian, "epsilon")
#'
#' ## Classification
#' pred <- bcnsmclass(fit_gaussian, wdbc[-id, 1:5], wdbc$Diagnosis[id])
#' table(pred = pred, true = wdbc$Diagnosis[-id])
#' }
#'
bcnsmclass <- function(object, test, cl){

  cl <- factor(cl)
  y <- as.matrix(test)
  d <- object$d
  pk <- prop.table(table(cl))

  if (object$association == "non-associative"){
    Gamma <- diag(d)
  }else{
    Gamma <- get(tolower(object$association))(d)$Gamma(round(object$gamma, 4))
  }

  mu <- object$fitted.values$mu
  sigma <- object$fitted.values$sigma
  lambda <- object$fitted.values$lambda
  nu <- object$fitted.values$nu

  out <- matrix(NA, nrow = nrow(y), ncol = length(pk))
  i <- 1
  for(k in levels(cl)){

    out[, i] <- dbcnsm(as.matrix(y), mu[which(cl == k), ][1, ], sigma, lambda, nu, Gamma,
                       copula = object$copula, delta = object$delta, margins = object$margins) * pk[i]

    i <- i + 1

  }

  aux <- function(x){
    which.max(x)
  }

  levels(cl)[apply(out, 1, which.max)]

}
