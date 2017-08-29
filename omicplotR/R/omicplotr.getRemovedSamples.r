#' Display removed samples.
#'
#' Generates table of removed samples (columns) by \code{omicplotr.filter}.
#'
#' @param x unfiltered data.
#' @param x.clr output from \code{omicplotr.clr}.
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
