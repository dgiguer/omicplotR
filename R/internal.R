#' Internal omicplotR functions.
#'
#' Internal functions for the app are not intended to be called by the
#'  user. Please refer to \link{omicplotr.run}.
#'
#' @author Daniel Giguere, Greg Gloor.
#'
#' @aliases omicplotr.anosim omicplotr.clr omicplotr.colvec omicplotr.filter
#'  omicplotr.colouredPCA omicplotr.getRemovedSamples
#'  omicplotr.getRemovedFeatures omicplotr.metadataFilter omicplotr.permanova
#'  omicplotr.report
#'
#' @export
#'
#' @examples
#' # These functions are not meant to be called by the user. However, some
#' # functions can be useful.
#'
#' data(otu_table)
#' # Filter out your OTU table to remove rows that sum to less than 10.
#'
#' d <- omicplotr.filter(otu_table, min.reads = 10)
#'
#' @import ALDEx2 igraph markdown propr
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

omicplotr.clr <- function(data, var.filt = 0) {

    #replace zeros
    #catch error if no 0 present
    if (any(data == 0)) {
        data.0 <- cmultRepl(t(data), label = 0, method = "CZM")
    } else {
        data.0 <- t(data)
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

omicplotr.colouredPCA <- function(data, colourvector, scale = 0,
    arrows = TRUE, taxonomy = NULL, show.taxonomy = FALSE,
    tax.level = 6, removenames = FALSE, names.cex = 1.0,
    points.cex = 0.8, main = "Principal Component Analysis Biplot")
    {
        x.var <- sum(data$sdev ^ 2)
        PC1 <- paste("PC 1 Variance: %",
        round(sum(data$sdev[1] ^ 2) / x.var * 100, 1))
        PC2 <- paste("PC 2 Variance: %",
        round(sum(data$sdev[2] ^ 2) / x.var*100, 1))

        points <- c(rep(".", length(dimnames(data$rotation)[[1]])))
        if ((isTRUE(show.taxonomy)) & (!is.null(taxonomy)))
        {
            genus <- sapply(strsplit(as.character(taxonomy),
            ";"), "[", tax.level)
        }

        col = c("black", rgb(0, 0, 0, 0.4))

        #remove sample names
        if (isTRUE(removenames)) {
            xlabs <- c(rep(".", length(dimnames(data$x)[[1]])))
            names.cex = 5*names.cex
            points.cex = points.cex
        } else {
            xlabs = unlist(dimnames(data$x)[1])
            names.cex = names.cex
            points.cex = points.cex
        }

        #if taxonomy is present, plot using taxonomies
        if (!isTRUE(show.taxonomy)) {
            points <- points
        } else {
            if (isTRUE(show.taxonomy)) {
                points <- genus
            } else {
                points <- points
            }
        }

        coloredBiplot(
            data,
            col = col,
            main = main,
            cex.main = 1.5,
            cex = c(names.cex, points.cex),
            xlabs.col = colourvector,
            var.axes = arrows,
            scale = scale,
            xlab = PC1,
            ylab = PC2,
            xlabs = xlabs,
            ylabs = points
        )
    }


    omicplotr.getRemovedFeatures <- function(x,
        x.clr) {

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

            #include OTUs as column because renderDataTable
            # doesn't include rownames
            output <- cbind(OTUs = rownames(output), output)

            return(output)

        }


        omicplotr.getRemovedSamples <- function(x, x.clr)
        {

            #get rid of taxonomy column if it exists
            x$taxonomy <- NULL

            #get original samples and retained samples
            names <- colnames(x)
            samples <- dimnames(x.clr$x)[[1]]

            missing <- c()

            #samples
            for (i in 1:length(names)) {
                if (!(names[i] %in% samples)) {
                    missing[i] <- names[i]
                }
            }

            missing <- missing[!is.na(missing)]

            if (is.null(missing)) {
                output <- missing
            } else {
                output <- x[, missing]

                #work around for including rownames in
                # dataTableOutput
                output <- cbind(OTUs = rownames(output),
                output)
            }

            return(output)

        }

        omicplotr.metadataFilter <- function(data, meta, column, values) {

            #order the metadata
            meta <- meta[order(rownames(meta)), ]

            #get rid of blank values
            s <- values[values != ""]

            #make list of samples
            meta.filt <- rownames(meta[which(meta[[column]]
                %in% s), ])

                #filter the data based on which samples remain
                # in metadata
                x <- data[which(colnames(data) %in% meta.filt)]

                if (is.null(data$taxonomy)) {
                    x <- x
                } else {
                    x$taxonomy <- data$taxonomy
                }

                return(x)
            }


            omicplotr.permanova <- function(x, conds) {

                d.clr <- apply(x, 1, function(x) {log(x) - mean(log(x))})

                dist.clr <- dist(d.clr)

                conds <- data.frame(conds)

                colnames(conds) <- "grp"

                test <- adonis(dist.clr~grp, data=conds, method="euclidean",
                    permutations=10000)

                    return(test$`Pr(>F)`)

                }



                omicplotr.report <- function(x.all, feature = NULL) {

                    features <-
                    rownames(x.all)[which(abs(x.all$effect) > 1 &
                    x.all$we.eBH
                    < 0.1)]

                    #effect sizes
                    effect <- x.all$effect[which(abs(x.all$effect)
                    > 1 & x.all$we.eBH < 0.1)]

                    #diff btw and diff win
                    diff.btw <-
                    x.all$diff.btw[which(abs(x.all$effect) > 1 &
                    x.all$we.eBH
                    < 0.1)]
                    diff.win <-
                    x.all$diff.win[which(abs(x.all$effect) > 1 &
                    x.all$we.eBH
                    < 0.1)]

                    #table of of significantly different OTUs with
                    # BH fdr < 0.1.
                    sig.output <- cbind(features, effect, diff.btw, diff.win)

                    }



                    omicplotr.filter <- function(data, min.reads = 0, min.count
                        = 0, min.sum = 0, min.prop = 0, max.prop = 1) {

                            if (is.null(data$taxonomy)) {
                                taxCheck <- TRUE
                            } else {
                                taxCheck <- FALSE
                            }

                            if (isTRUE(taxCheck)) {
                                #order for colouring
                                x <- data[order(colnames(data))]

                                #filter by min reads per sample
                                data.0 <- x[,which(apply(x,2,sum) > min.reads)]

                                #filter by min count per feature
                                data.1 <-  data.0[which(apply(data.0, 1, max) >=
                                min.count), ]

                                #filter by sum of count per feature
                                data.2 <- data.1[which(apply(data.1, 1, sum) >
                                min.sum), ]

                                #filter by proportion
                                d.frac <- apply(data.2, 2,
                                    function(x){x/sum(x)})
                                x.filt <- data.2[ (which((apply(d.frac,1,max) >
                                min.prop) &
                                (apply(d.frac,1,max) < max.prop))),]

                            } else {
                                #order for colouring
                                tax <- data$taxonomy
                                x <- data[order(colnames(data[1:(ncol(data) -
                                 1)]))]
                                x$taxonomy <- tax

                                #filter by min reads per sample
                                data.0 <-  data.frame(x[,which(apply
                                    (x[, 1:(ncol(x) - 1)], 2, sum) >=
                                     min.reads)], check.names = FALSE)

                                    #put tax back in
                                    data.0$taxonomy <- tax
                                    #filter by min count per feature
                                    data.1 <-  data.frame(data.0[which(apply
                                        (data.0[,1:(ncol(data.0) - 1)], 1, max)
                                         >= min.count),], check.names = FALSE)

                                        #filter by sum of count per feature
                                        data.2 <- data.frame(data.1[which(apply
                                            (data.1[, 1:(ncol(data.1) - 1)], 1,
                                             max) >= min.sum), ], check.names =
                                              FALSE)

                                            #filter by proportion
                                            d.frac <-
                                             apply(data.2[,1:(ncol(data.2) -
                                            1)], 2,
                                            function(x){x/sum(x)})
                                            x.filt <- data.2[
                                             (which((apply(d.frac,1,max)
                                            > min.prop) &
                                            (apply(d.frac,1,max) < max.prop))),]
                                        }
                                        return(x.filt)
                                    }


  omicplotr.colvec <- function(data, meta, column, type = 3) {

    df <- as.data.frame(unlist(dimnames(data$x)[1]))

    colnames(df) <- "samples"
    #add empty column (which is gonna be the ample list for
    # colouring)
    df$correct <- NA
    #introduces NA where metadata is not present for sample
    #the important thing is that this keeps the order of the data
    for (i in 1:length(df$samples)) {
      if (df$samples[i] %in% rownames(meta)) {
        df$correct[i] <- as.character(df$samples[i])
      } else {
        df$correct[i] <- NA
      }
    }

    #colour column
    df$colour <- NA

    #if correct != present, colour black
    #this means metadata was missing for data point
    for (i in 1:length(df$correct)) {
      if (is.na(df$correct[i])) {
        df$colour[i] <- "black"
      }
    }

    #get unique features of columns in metadata
    u <- unique(meta[[column]])
    #get rid of NAs
    u <- u[!is.na(u)]

    #default colour palette
    colours <- c("indianred1", "steelblue3",  "skyblue1",
    "mediumorchid", "royalblue4", "olivedrab3", "pink", "#FFED6F",
    "mediumorchid3", "ivory2"
  )

  #colour the colour column for up to ten colours
  for (i in 1:length(df$correct)) {
    if (is.na(meta[column][df$correct[i], ])) {
      df$colour[i] <- "black"
    } else if (meta[column][df$correct[i], ] ==
      as.character(u[1])) {
        df$colour[i] <- colours[1]
      } else if (meta[column][df$correct[i], ] ==
        as.character(u[2])) {
          df$colour[i] <- colours[2]
        } else if (meta[column][df$correct[i], ] ==
          as.character(u[3])) {
            df$colour[i] <- colours[3]
          } else if (meta[column][df$correct[i], ] ==
            as.character(u[4])) {
              df$colour[i] <- colours[4]
            } else if (meta[column][df$correct[i], ] ==
              as.character(u[5])) {
                df$colour[i] <- colours[5]
              } else if (meta[column][df$correct[i], ] ==
                as.character(u[6])) {
                  df$colour[i] <- colours[6]
                } else if (meta[column][df$correct[i], ] ==
                  as.character(u[7])) {
                    df$colour[i] <- colours[7]
                  } else if (meta[column][df$correct[i], ] ==
                    as.character(u[8])) {
                      df$colour[i] <- colours[8]
                    } else if (meta[column][df$correct[i], ] ==
                      as.character(u[9])) {
                        df$colour[i] <- colours[9]
                      } else if (meta[column][df$correct[i], ] ==
                        as.character(u[10])) {
                          df$colour[i] <- colours[10]
                        }
                      }

                      #different colouringtype
                      if (type == 1) {

                        #TODO get rid of this if using command line
                        validate(need(is.integer(meta[[column]]), "Nonnumeric
                         data detected. Please click 'Nonnumeric metadata'"))

                        #get rid of old column made
                        df$colour <- NULL

                        #make age column
                        df$age <- NA

                        #get the age of each sample
                        for (i in 1:length(df$correct)) {
                          if (is.na(meta[column][df$correct[i], ])) {
                            df$age[i] <- NA
                          } else {
                            df$age[i] <- meta[column][df$correct[i], ]
                          }
                        }

                        #order data frame in order
                        df.sort <- df[order(df$age), ]

                        df.sort$col <- NA

                        #colour each unique number seperately
                        df.u <- unique(df.sort$age)

                        len <- length(df.u)

                        #generate colour function
                        c <- colorRampPalette(c("red", "blue"))(len)

                        #colour each unique age
                        for (i in 1:length(df.sort$age)) {
                          df.sort$col[df.sort$age == df.u[i]] <- c[i]
                        }

                        #colour the NAs black. these are missing data points.
                        for (i in 1:length(df.sort$age)) {
                          if (is.na(df.sort$col[i])) {
                            df.sort$col[i] <- "#000000"
                          }
                        }

                        #reorder the df by sample now so it matches the data
                        df.sort <- df.sort[order(df.sort$samples), ]

                        #voila, c'est le colour!
                        sorted <- c()
                        sorted$col <- df.sort$col

                      } else if (type == 2) {

                        validate(need(is.integer(meta[[column]]), "Nonnumeric
                         data detected. Please click 'Nonnumeric metadata'"))

                        #quartile
                        df$colour <- NULL

                        #make age column
                        df$age <- NA

                        #get the age of each sample
                        for (i in 1:length(df$correct)) {
                          if (is.na(meta[column][df$correct[i], ])) {
                            df$age[i] <- NA
                          } else {
                            df$age[i] <- meta[column][df$correct[i], ]
                          }
                        }

                        df$col <- NA

                        c <- colorRampPalette(c("green", "blue"))(4)

                        #make quartile rank column
                        sorted <-
                        within(df, col <-
                          as.integer(cut(
                            df$age, quantile(df$age, na.rm = TRUE),
                            include.lowest = TRUE
                          )))

                          #replace quartile rank with colours
                          sorted$col[sorted$col == 1] <- c[1]
                          sorted$col[sorted$col == 2] <- c[2]
                          sorted$col[sorted$col == 3] <- c[3]
                          sorted$col[sorted$col == 4] <- c[4]

                          #replace any NAs with black
                          for (i in 1:length(sorted$col)) {
                            if (is.na(sorted$col[i])) {
                              sorted$col[i] <- "#000000"
                            }
                          }

                        } else if (type == 3) {
                          sorted <- c()
                          sorted$col <- df$colour
                        } else {
                          return("")
                        }
                        return(sorted$col)
                      }
