#' Calculate \code{anosim} between conditions.
#'
#' Wraps \code{anosim} from \code{vegan} package to calculate anosim
#'  p-value between conditions from cenetered-log ratio values between two
#'  conditions.
#'
#' @details
#' Any zeros are replaced using the CZM method using \code{cmultRepl} from
#'  \code{zCompositions}. \code{dist} computes the distance matrix from the
#'  cenetered-log ratio.
#'
#' @param x Filtered data output from \code{omicplotr.filter}. Columns must
#'  in ordered of condition being tested (i.e., first 4 columns condition x and
#'  next 4 columns condition 4).
#' @param conds A character vector of group labels. See example.
#'
#' @return Returns significance from permutation.
#'
#' @seealso \code{\link[vegan]{anosim}},
#'  \code{\link[zCompositions]{cmultRepl}}, \code{\link[ALDEx2]{selex}}
#'
#' @examples
#' # For selex example data set where first 7 samples are non-selected, the last
#' # 7 are selected.
#' # conds <- c(rep("NS", 7), rep("S", 7))
#' # p <- omicplotr.anosim(selex, conds)
#'
#' @export
#'

omicplotr.anosim <- function(x, conds) {
    #get rid of tax column
    x$taxonomy <- NULL

    #replace zeros
    d <- cmultRepl(t(x), label = 0, method = "CZM")

    #calculate clr
    d.clr <- apply(d, 1, function(x) {
        log(x) - mean(log(x))
      })

    #dissimilarity matrix
    d.dist <- dist(d.clr)

    #calculate anosim, requires
    anosim(d.dist, conds, permutations = 999)

    #p-value of between conditions
    return(anosim$signif)
}
