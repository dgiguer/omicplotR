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
#' @examples
#' # To run omicplotR, use `omicplotr.run()` function. This will open the shiny
#' # app in your default browser.
#'
#' # omicplotr.run()
#'
#' # After running omicplotr.run(), you can filter your data and make a PCA
#' # biplot. Internal functions are not intended to be called by the user.
#'
#' data(otu_table)
#' d.filt <- omicplotr.filter(otu_table)
#'
#' @value
#' Will launch the Shiny App. There is no value.
#'
#' @import ALDEx2 igraph markdown propr
#'

omicplotr.run <- function() {
    
    appDir <- system.file("shiny-app", package = "omicplotR")
    if (appDir == "") {
        stop("Could not find directory. Try reinstalling omicplotR", call = FALSE)
    }
    
    shiny::runApp(appDir)
}
