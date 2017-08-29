#' Calculate \code{anosim} between conditions.
#'
#' Wrapper for \code{vegan} package to calculate anosim p-value between
#'  conditions from cenetered-log ratio values between two conditions.  
#'
#' @param x Filtered data output from \code{omicplotr.filter}. Columns must
#'  in ordered of condition being tested (i.e., first 4 columns condition x and
#'  next 4 columns condition 4).
#' @param conds a character vector of group labels.
#'
#' @examples
#' For dataset where first 7 columns are condition x, next 7 are condition:
#' conds <- c(rep("group1", 7). rep("group2", 7))
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
    ansoim(d.dist, conds, permutations = 999)

    #p-value of between conditions
    return(anosim$signif)
}
