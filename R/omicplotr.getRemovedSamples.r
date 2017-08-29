#' Display removed samples.
#'
#' Generates table of removed samples (columns) by \code{omicplotr.filter}.
#'
#' @param x Unfiltered data.
#' @param x.clr Output from \code{omicplotr.clr}.
#'
#' @seealso \code{\link{omicplotr.clr}}
#'
#' @author Daniel Giguere
#'
#' @export

omicplotr.getRemovedSamples <- function(x, x.clr) {

    #get rid of taxonomy column if it exists
    x$taxonomy <- NULL

    #get original samples and retained samples
    names <- colnames(x)
    samples <- dimnames(x.clr$x)[[1]]

    missing <- c()

    #samples
    for (i in 1:length(names)) {
        if (!(names[i] %in% samples)) {
            missing[i] <- names[i]
        }
    }

    missing <- missing[!is.na(missing)]

    if (is.null(missing)) {
        output <- missing
        } else {
            output <- x[, missing]

            #work around for including rownames in dataTableOutput
            output <- cbind(OTUs = rownames(output), output)
        }

        return(output)

    }
