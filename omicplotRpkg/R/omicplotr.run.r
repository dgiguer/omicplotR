# in progress
#' Run \code{omicplotR} \code{Shiny} app.
#'
#' Launches \code{Shiny} for \code{omicplotR} in default browser.
#'
#' @param x Filtered data output from \code{omicplotr.filter}. Columns must
#'  in ordered of condition being tested (i.e., first 4 columns condition x and
#'  next 4 columns condition 4).
#' @param conds a character vector of group labels.
#'
#' @examples
#' For dataset where first 7 columns are condition x, next 7 are condition:
#' conds <- c(rep("group1", 7). rep("group2", 7))
#' @export
#'

omicplotr.run <- function() {

    appDir <- system.file("inst", package = "omicplotRpkg")
    if (appDir == "") {
        stop("Could not find directory. Try again", call=FALSE)
    }

    shiny::runApp(appDir)
}
