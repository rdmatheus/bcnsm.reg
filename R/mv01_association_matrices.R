## Nonassociative ----------------------------------------------------------------------------------
nonassociative <- function(d) {
  ans <- list()

  ans$npar <- 0L

  ans$start <- function(x) NULL

  ans$Gamma <- function(gamma = NULL, ...) diag(d)

  ans$lower <- NULL
  ans$upper <- NULL

  ans$name <- "non-associative"

  class(ans) <- "association"
  ans
}


## Unstructured ------------------------------------------------------------------------------------
unstructured <- function(d) {
  ans <- list()

  ans$npar <- d * (d - 1) / 2

  ans$start <- function(x){

    G <- sin(0.5 * pi * stats::cor(x, method = "kendall"))
    G[lower.tri(G, diag = TRUE)] <- NA

    gamma <- as.numeric(stats::na.exclude(c(t(G))))

  }

  ans$Gamma <- function(gamma) {
    class(gamma) <- "dist"
    attr(gamma, "Size") <- d
    gamma <- as.matrix(gamma)
    diag(gamma) <- 1
    gamma
  }

  EPS <- sqrt(.Machine$double.eps)
  ans$lower <- rep(-1 + EPS, d * (d - 1) / 2)
  ans$upper <- rep(1 - EPS, d * (d - 1) / 2)

  ans$name <- "unstructured"

  class(ans) <- "association"
  ans
}

### Uniform (exchangeable) -------------------------------------------------------------------------
uniform <- function(d) {
  ans <- list()

  ans$npar <- 1

  ans$start <- function(x) stats::median(sin(0.5 * pi * stats::cor(x, method = "kendall"))[lower.tri(diag(d))])

  ans$Gamma <- function(gamma) {
    out <- matrix(gamma, d, d)
    diag(out) <- 1
    out
  }

  EPS <- sqrt(.Machine$double.eps)
  ans$lower <- -1 + EPS
  ans$upper <- 1 - EPS

  ans$name <- "uniform"

  class(ans) <- "association"
  ans
}



