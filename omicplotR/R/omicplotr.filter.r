#' Filter data frame by read counts.
#'
#' \code{omicplotr.filter} removes columns and/or rows of table below threshold.
#' This must be used after \code{omicplotr.metadataFilter}.
#'
#' @param x counts table as data frame.
#' @param min.reads removes any columns with sum less than or equal to input.
#' @param min.count removes any rows with maximum lower than input.
#' @param min.sum removes any rows with sum lower than input.
#' @param min.prop minimum proportional abundance in any sample.
#' @param max.prop maximum proportional abundance in any sample.
#'
#' @examples
#' # The 'reads' data frame must have unique row and column names, and may
#' # include a taxonomy column as the last column, and look like this
#' # (from \code{ALDEx2} documentation):
#' #
#' #              T1a T1b  T2  T3  N1  N2  Nx taxonomy
#' #   Gene_00001   0   0   2   0   0   1   0 Bacteria;Bacteroidetes;...
#' #   Gene_00002  20   8  12   5  19  26  14 Bacteria;Firmicutes;...
#' #   Gene_00003   3   0   2   0   0   0   1 Bacteria;Proteobacteria;...
#' #   Gene_00004  75  84 241 149 271 257 188 Bacteria;Bacteroidetes;...
#' #   Gene_00005  10  16   4   0   4  10  10 Bacteria;unclassified;...
#' #   Gene_00006 129 126 451 223 243 149 209 Bacteria;Firmicutes;...
#' #       ... many more rows ...
#' #
#' #   The taxonomy column is not required.
#' #
#' @export

omicplotr.filter <- function(x, min.reads = 0, min.count = 0, min.sum = 0,
     min.prop = 0, max.prop = 1) {

        if (is.null(x$taxonomy)) {
            taxCheck <- TRUE
            } else {
            taxCheck <- FALSE
            }

            if (isTRUE(taxCheck)) {
                #order for colouring
                x <- x[order(colnames(x))]

                #filter by min reads per sample
                data.0 <- x[,which(apply(x,2,sum) > min.reads)]

                #filter by min count per feature
                data.1 <-  data.0[which(apply(data.0, 1, max) >= min.count), ]

                #filter by sum of count per feature
                data.2 <- data.1[which(apply(data.1, 1, sum) > min.sum), ]

                #filter by proportion
                d.frac <- apply(data.2, 2, function(x){x/sum(x)})
                x.filt <- data.2[ (which((apply(d.frac,1,max) > min.prop) &
                (apply(d.frac,1,max) < max.prop))),]

                } else {
                    #order for colouring
                    tax <- x$taxonomy
                    x <- x[order(colnames(x[1:(ncol(x) - 1)]))]
                    x$taxonomy <- tax

                    #filter by min reads per sample
                    data.0 <-  data.frame(x[,which(apply
                        (x[, 1:(ncol(x) - 1)], 2, sum) >= min.reads)],
                        check.names = FALSE)

                    #put tax back in
                    data.0$taxonomy <- tax
                    #filter by min count per feature
                    data.1 <-  data.frame(data.0[which(apply
                    (data.0[,1:(ncol(data.0) - 1)], 1, max) >= min.count),],
                    check.names = FALSE)

                    #filter by sum of count per feature
                    data.2 <- data.frame(data.1[which(apply
                    (data.1[, 1:(ncol(data.1) - 1)], 1, max) >= min.sum), ],
                    check.names = FALSE)

                    #filter by proportion
                    d.frac <- apply(data.2[,1:(ncol(data.2) - 1)], 2,
                    function(x){x/sum(x)})
                    x.filt <- data.2[ (which((apply(d.frac,1,max) > min.prop) &
                    (apply(d.frac,1,max) < max.prop))),]
                }
                return(x.filt)
            }
