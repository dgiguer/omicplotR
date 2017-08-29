#' Generate summary report of significantly different features.
#'
#' Displays table of features found to be significantly different by
#'  \code{ALDEx2}, when effect size > 1 and expected Benjamini Hochberg
#'  corrected p-values. In progress and will show more.
#'
#' @param x.all output from \code{ALDEx2}
#' @param feature optional; report for single feature/OTU
#'
#' @examples
#'
#' @export
#'
omicplotr.report <- function(x.all, feature = NULL) {

    features <- rownames(x.all)[which(abs(x.all$effect) > 1 & x.all$we.eBH
    < 0.1)]

    #effect sizes
    effect <- x.all$effect[which(abs(x.all$effect) > 1 & x.all$we.eBH < 0.1)]

    #diff btw and diff win
    diff.btw <- x.all$diff.btw[which(abs(x.all$effect) > 1 & x.all$we.eBH
    < 0.1)]
    diff.win <- x.all$diff.win[which(abs(x.all$effect) > 1 & x.all$we.eBH
    < 0.1)]

    #table of of significantly different OTUs with BH fdr < 0.1.
    sig.output <- cbind(features, effect, diff.btw, diff.win)

}
