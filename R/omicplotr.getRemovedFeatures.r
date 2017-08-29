#' Display removed features.
#'
#' Generates table of features (rows) removed by \code{omicplotr.filter}.
#'
#' @param x Unfiltered data.
#' @param x.clr Output from \code{omicplotr.clr}.
#'
#' @seealso \code{\link{omicplotr.clr}}
#'
#' @author Daniel Giguere
#'
#' @export

omicplotr.getRemovedFeatures <- function(x, x.clr) {

    #get rid of taxonomy column if it exists.
    x$taxonomy <- NULL

    #get original samples and retained samples
    otu.names <- rownames(x)
    otus <- dimnames(x.clr$rotation)[[1]]

    otu.missing <- c()

    #samples
    for (i in 1:length(otu.names)) {
        if (!(otu.names[i] %in% otus)) {
            otu.missing[i] <- otu.names[i]
        }
    }

    otu.missing <- otu.missing[!is.na(otu.missing)]

    output <- x[otu.missing,]

    #include OTUs as column because renderDataTable doesn't include rownames
    output <- cbind(OTUs = rownames(output), output)

    return(output)

    }
