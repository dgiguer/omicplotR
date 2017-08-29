#' Calculate PERMANOVA between conditions.
#'
#' Wrapper for \code{vegan} package to calculate PERMANOVA between
#'  conditions from cenetered-log ratio values with \code{adonis}.
#'
#' @param x Filtered data output from \code{omicplotr.filter}. columns must
#'  in ordered of condition being tested (i.e., first 7 columns condition x and
#'  next 7 columns condition 2).
#' @param conds A character vector of group labels. See example.
#'
#' @examples
#' # For selex dataset.
#' conds <- c(rep("NS", 7), rep("S", 7))
#'
#' @seealso \code{\link[vegan]{adonis}}, \code{\link[ALDEx2]{selex}}
#'
#' @export
#'

omicplotr.permanova <- function(x, conds) {

    d.clr <- apply(x, 1, function(x) {log(x) - mean(log(x))})

    dist.clr <- dist(d.clr)

    conds <- data.frame(conds)

    colnames(conds) <- "grp"

    test <- adonis(dist.clr~grp, data=conds, method="euclidean",
    permutations=10000)

    return(test$`Pr(>F)`)

}
