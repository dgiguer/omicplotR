#' Calculate CLR and generate prcomp object for compositional biplot
#'
#' \code{omciplotr.clr} replaces zeros using CZM method using \code{cmultRepl}
#'  from \code{zCompositions}, calculates CLR, filters by variance and
#'  generates a \code{prcomp} object from data.
#'
#' @details If filtering by with \code{omicplotr.filter} or
#'  \code{omicplotr.metadataFilter}, it must done before this step.
#'
#' @param x filtered data from \code{omicplotr.filter}.
#' @param var.filt variance cutoff to be applied after clr is computed,
#'  default, 0.
#'
#' @return Returns an object of class "prcomp".
#'
#' @seealso \code{\link[stats]{prcomp}}.
#'
#' @export
#'

# replace zeros, calculate CLR and convert to prcomp object
omicplotr.clr <- function(x, var.filt = 0) {

    #replace zeros
    #catch error if no 0 present
    if (any(x == 0)) {
        data.0 <- cmultRepl(t(x), label = 0, method = "CZM")
        } else {
            data.0 <- t(x)
        }

        #calculate CLR
        x.clr <-
        t(apply(data.0, 1, function(x) {
            log(x) - mean(log(x))
            }))

            #filter by variance
            if (var.filt == 0) {
                x.clr.var <- x.clr
                } else {
                    var.clr <- apply(x.clr, 2, var)
                    names.hvar <- names(var.clr)[which(var.clr > var.filt)]
                    x.clr.var <- x.clr[,names.hvar]
                }

                #generate prcomp object
                x.pcx <- prcomp(x.clr.var)
            }
