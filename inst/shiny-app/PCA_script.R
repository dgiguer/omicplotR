#script for making reproducing the PCA plot. values chosen for filtering in omicplotR are listed above
#subsetting by metadata is intended for coloured biplots so it isn't included here
#all omicplotr functions can be found in inst/shiny-app/internal.r

library(zCompositions)
library(ALDEx2)
library(omicplotR)

##################################################################
#                                                                #
# You must change the file variable to your ownfilename           #
#                                                                #
#                                                                #
##################################################################
#change filename.txt to your filename, make sure it follow input requirements
file <- "filename.txt"

#read file
data <- read.table(
  file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  row.names = 1,
  check.names = FALSE,
  comment.char = "",
  na.strings = "")

#filter data
x.filt <- omicplotr.filter(data, min.reads = min.reads, min.count = min.count, min.prop = min.prop, max.prop = max.prop, min.sum = min.sum)

#check for taxonomy column
if (is.null(data$taxonomy)) {
  taxCheck <- TRUE
} else {
  taxCheck <- FALSE
}

if (isTRUE(taxCheck)) {
  x.filt <- x.filt
} else {
  x.filt$taxonomy <- NULL
}

#calculates clr and generates prcomp object
data.pr <- omicplotr.clr(x.filt, var.filt)

#calculate variance for x/y axis of plot
x.var <- sum(data.pr$sdev ^ 2)
PC1 <- paste("PC 1 Variance: %", round(sum(data.pr$sdev[1] ^ 2) / x.var * 100, 1))
PC2 <- paste("PC 2 Variance: %", round(sum(data.pr$sdev[2] ^ 2) / x.var*100, 1))

##############################################################################

#checks for biplot

#take taxonomy column from data if it exists
taxonomy <- x.filt$taxonomy

#if Show taxonomy checkbox is clicked, take chosen level or genus by default
if (isTRUE(taxoncheck)) {
    taxselect <- as.numeric(input$taxlevel)
} else {
    taxselect <- 6
}

#get genus (or other level)
genus <- vapply(strsplit(as.character(taxonomy), "[[:punct:]]"), "[", taxselect, FUN.VALUE=character(1))

#do checks for arrows
if (isTRUE(arrowcheck)) {
  arrows = FALSE
} else {
  arrows = TRUE
}
    #if taxonomy is is null, use periods as points for features.
    #if taxoncheckbox is checked, show genus instead of points
  if (isTRUE(taxoncheck)) {
    points <- genus
  } else {
    points <- c(rep(".", length(dimnames(data.pr$rotation)[[1]])))
  }

#remove sample names if checkbox clicked
if (isTRUE(removenames)) {
  xlabs <- c(rep(".", length(dimnames(data.pr$x)[[1]])))
  size = c(5, 0.8)
} else {
  xlabs = unlist(dimnames(data.pr$x)[1])
  size = c(1.0, 0.8)
}

##############################################################################

#colour of points (features black, samples red)
col = c("black", rgb(0, 0, 0, 0.2))

#add or change options for biplot.
biplot(
  data.pr,
  main = "Principal Component Analysis",
  cex.main = 1.5,
  cex = size,
  col = col,
  scale = scale.slider,
  var.axes = arrows,
  xlab = PC1,
  ylab = PC2,
  xlabs = xlabs,
  ylabs = points
)
