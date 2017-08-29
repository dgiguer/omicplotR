#' Calculate PERMANOVA between conditions.
#'
#' Wrapper for \code{vegan} package to calculate PERMANOVA between
#'  conditions from cenetered-log ratio values.
#'
#' @param x filtered data output from \code{omicplotr.filter}. columns must
#'  in ordered of condition being tested (i.e., first 4 columns condition x and
#'  next 4 columns condition 4).
#' @param conds a character vector of group labels.
#'
#' @examples
#' for dataset where first 7 columns are condition x, next 7 are condition y.
#' conds <- c(rep("group1", 7). rep("group2", 7))
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
