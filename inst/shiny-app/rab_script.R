#script for reproducing the dendrogram and relative abundance stacked barplots. values chosen for filtering in omicplotR are listed above
#all omicplotr functions can be found in inst/shiny-app/internal.r

library(zCompositions)
library(ALDEx2)
library(omicplotR)

##################################################################
#                                                                #
# You must change the file variable to your own file name        #
#                                                                #
#                                                                #
##################################################################

#change filename.txt to your filename, make sure it follow input requirements
file <- "filename.txt"

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

#check for taxonomy column, get taxonomy if present
if (is.null(x.filt$taxonomy)) {
  print("Taxonomy column required")
} else {
  d <- x.filt[, seq_along(x.filt) - 1]
  taxon <- x.filt[(dim(x.filt)[2])]
}

# Get genera
genera <- c()
for(i in (seq_len(nrow(taxon)))) {
  genera <- c(genera, vapply(strsplit(as.character(taxon[i, ]), "[[:punct:]]"),
  "[", 6, FUN.VALUE=character(1)))
}

# sum counts by name
d.agg <- aggregate(d, by=list(genera), FUN=sum)
tax.agg <- d.agg$Group.1
d.agg$Group.1 <- NULL

# convert to abundances
d.prop <- apply(d.agg, 2, function(x){x/sum(x)})

#filters by abundance (slider bar)
d.abund <- d.agg[apply(d.prop, 1, max) > abund,]
tax.abund.u <- tax.agg[apply(d.prop, 1, max) > abund]

d.abund <- t(cmultRepl(t(d.abund), label=0, method="CZM"))

# get proportions of the filtered data for plotting below
# in log-ratio speak, you are re-closing your dataset
d.P.u <- apply(d.abund, 2, function(x){x/sum(x)})

# order by OTU abundances
new.order <- rownames(d.P.u)[order(apply(d.P.u, 1, sum), decreasing=TRUE)]
tax.abund <- tax.abund.u[order(apply(d.P.u, 1, sum), decreasing=TRUE)]
d.P <- d.P.u[new.order, ]
d.clr <- apply(d.P, 2, function(x){log2(x) - mean(log2(x))})

#distance matrix
method.m <- c("euclidean", "maximum", "manhattan")
method.mchoice <- method.m[dist]

dist.d.clr <- dist(t(d.clr), method=method.mchoice)

#get clustering method for dendrogram
method <- c("complete", "single", "ward.D2")
method.choice <- method[clust]

clust.d <- hclust(dist.d.clr, method=method.choice)

#plot the dendrogram (option for wide pdfs below)
plot(as.dendrogram(clust.d), main=NULL, cex=0.8, xlab="", xpd = TRUE)

#change colours if you'd like
colours <- c("steelblue3", "skyblue1", "indianred", "mediumpurple1",
"olivedrab3", "pink", "#FFED6F", "mediumorchid3",
"ivory2", "tan1", "aquamarine3", "#C0C0C0", "royalblue4",
"mediumvioletred", "#999933", "#666699", "#CC9933", "#006666",
"#3399FF", "#993300", "#CCCC99", "#666666", "#FFCC66",
"#6699CC", "#663366", "#9999CC", "#CCCCCC", "#669999",
"#CCCC66", "#CC6600", "#9999FF", "#0066CC", "#99CCCC",
"#999999", "#FFCC00", "#009999", "#FF9900", "#999966",
"#66CCCC", "#339966", "#CCCC33", "#EDEDED")

#relative abundance barplot
plot.new()

#change sizing so barplot and legend fit
par(fig=c(0,1,0,1), new = TRUE)
par(fig=c(0,0.80,0, 1), new = TRUE)

#stacked barplot
barplot(d.P[,clust.d$order], names.arg = clust.d$labels, space=0, col=colours, las=2, axisnames=TRUE, border = NA, xpd = TRUE)

#change location of legend
par(fig=c(0.8,1, 0, 1), new=TRUE)
par(xpd = TRUE)

#get legend colours
leg.col <- rev(colours[seq_len(nrow(d.P))])

#add legend
legend(x="center", legend=rev(tax.abund), col=leg.col, lwd=5, cex=0.8,
border=NULL)

######## option for wide pdfs

#approximate width for clear sample names, change if too wide/narrow.
#remove hash tags and enter into R.

######## dendrogram
#wid <- round(length(clust.d$labels)/5)
#pdf("omicplotR_dendrogram.pdf", width = wid)
#plot(as.dendrogram(clust.d), main=NULL, cex=0.8, xlab="", xpd = T)
#dev.off()

######## barplot
# wid <- round(length(clust.d$labels)/5)
# pdf("omicplotR_rab_barplot.pdf", width = wid)
# #change sizing so barplot and legend are on same graph
# par(fig=c(0,1,0,1), new = TRUE)
# par(fig=c(0,0.80,0, 1), new = TRUE)
#
# #stacked barplot
# barplot(d.P[,clust.d$order], names.arg = clust.d$labels, space=0, col=colours, las=2, axisnames=T, border = NA, xpd = T)
#
# #change location of legend
# par(fig=c(0.8,1, 0, 1), new=TRUE)
# par(xpd = TRUE)
#
# #get legend colours
# leg.col <- rev(colours[seq_len(nrow(d.P))])
#
# #add legend
# legend(x="center", legend=rev(tax.abund), col=leg.col, lwd=5, cex=0.8,
# border=NULL)
# dev.off()

########
################################################################################
