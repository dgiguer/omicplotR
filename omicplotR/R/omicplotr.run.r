# in progress
#' Run \code{omicplotR} \code{Shiny} app.
#'
#' Launches \code{Shiny} for \code{omicplotR} in default browser.
#' 
#' @export
#'

omicplotr.run <- function() {

    appDir <- system.file("shiny-app", package = "omicplotR")
    if (appDir == "") {
        stop("Could not find directory. Try reinstalling omicplotR", call=FALSE)
    }

    shiny::runApp(appDir)
}
