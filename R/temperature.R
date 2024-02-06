#' Monthly Climatological Data at Congonhas Airport
#'
#' Climatological data measured monthly at Congonhas Airport (longitude -46.6553, latitude -23.6267),
#'     Brazil, from January 2001 to December 2022. The data is collected by the NASA Prediction Of
#'     Worldwide Energy Resources (POWER) project and is also available through the R package
#'     \code{nasapower} (Sparks, 2018).
#'
#'
#' @format ## `temperature`
#'
#' A data frame with 264 rows and 7 columns:
#' \describe{
#'   \item{Year}{Year of measurement.}
#'   \item{Month}{Month of measurement.}
#'   \item{Temperature}{Monthly average air temperature at 2 meters above the surface of the earth,
#'       measured in degrees Celsius (\emph{°C}).}
#'   \item{Wind}{Average wind speed in the previous month, in \emph{m/s},}
#'   \item{Prec}{Bias-corrected average of total precipitation at the surface of the earth in the
#'       previous month, in \emph{mm/day}}
#'   \item{Time}{A POSIXct variable with the measurement date.}
#'   \item{Season}{Season of the year in which the measurement was taken.}
#' }
#'
#' @source The NASA Prediction Of Worldwide Energy Resources (POWER) project
#'     \url{}.
#'
#' @references
#'
#' Sparks, A. H. (2018). nasapower: a NASA POWER global meteorology, surface solar energy
#'     and climatology data client for R. \emph{The Journal of Open Source Software}, 3(30), 1035
#'
#'
"temperature"
