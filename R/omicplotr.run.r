#' Run \code{omicplotR} \code{Shiny} app.
#'
#' Launches \code{Shiny} for \code{omicplotR} in default browser. This
#'  graphical user interface facilitates visual exploration of datasets for
#'  both new and experienced R users alike.
#'
#' @author Daniel Giguere
#'
#' @export
#'
#' @import ALDEx2 igraph markdown propr
#'

omicplotr.run <- function() {

    appDir <- system.file("shiny-app", package = "omicplotR")
    if (appDir == "") {
        stop("Could not find directory. Try reinstalling omicplotR", call=FALSE)
    }

    shiny::runApp(appDir)
}
