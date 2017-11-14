#' Filter data by selected metadata.
#'
#' Filters samples in data to remain using values from associated metadata.
#'
#' @details This filter must be applied before filtering by count using
#'  \code{omicplotr.filter} or generating \code{prcomp} object with
#'  \code{omicplotr.clr}.
#'
#' @param data Unfiltered data.
#' @param meta Unfiltered metadata.
#' @param column Column name from metadata. Must exactly match.
#' @param values List of values to keep from metadata column. Must exactly match
#'  the values in the metadata (no extra characters - including spaces).
#'
#' @author Daniel Giguere
#'
#' @seealso \code{\link{omicplotr.filter}}, \code{\link{omicplotr.clr}}.
#'
#' @export

omicplotr.metadataFilter <- function(data, meta, column, values) {

    #order the metadata
    meta <- meta[order(rownames(meta)), ]

    #get rid of blank values
    s <- values[values != ""]

    #make list of samples
    meta.filt <- rownames(meta[which(meta[[column]] %in% s), ])

    #filter the data based on which samples remain in metadata
    x <- data[which(colnames(data) %in% meta.filt)]

    if (is.null(data$taxonomy)) {
      x <- x
    } else {
      x$taxonomy <- data$taxonomy
    }

    return(x)
}
