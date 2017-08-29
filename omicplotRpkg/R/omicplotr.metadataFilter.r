#' Filter data by selected metadata.
#'
#' Filters data by choosing values to remain using associated metadata.
#'
#' @param data unfiltered data.
#' @param meta unfiltered metadata.
#' @param column column name from metadata.
#' @param values list of values to keep from metadata column.
#'
#' @examples
#' values <- c("y", "n")
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

    return(x)
}
